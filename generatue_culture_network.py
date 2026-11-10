import numpy as np
import matplotlib.pyplot as plt
import random
import sys
import networkx as nx
from matplotlib.patches import Patch

seed = int(sys.argv[3])

np.random.seed(seed)
random.seed(seed)

# Define constants
circle_diameter = float(sys.argv[2])  # in mm
circle_radius = circle_diameter / 2
num_neurons = int(sys.argv[1])
neuron_diameter = 15e-3  # in mm
neuron_radius = neuron_diameter / 2
alpha = 0.1  # connection probability when axon intersects a dendritic tree

margin = 0.001

# Compute area of the circular culture
culture_area = np.pi * (circle_radius ** 2)
neuron_density = num_neurons / culture_area

# Dendritic tree parameters (in mm)
mu = 300e-3     # mean diameter in mm
sigma = 0.04    # std dev in mm

# Axon parameters (in mm)
sigma_l = 0.9
mean_length = 1.1
segment_length = 0.01
sigma_theta = 0.1

def is_valid_position(x, y, existing_positions):
    if x**2 + y**2 > circle_radius**2:
        return False
    for (ex, ey) in existing_positions:
        if np.sqrt((x - ex)**2 + (y - ey)**2) < 2 * neuron_radius:
            return False
    return True

def reflect_direction(direction, normal):
    direction = np.array(direction) / np.linalg.norm(direction)
    normal = np.array(normal) / np.linalg.norm(normal)
    return direction - 2 * np.dot(direction, normal) * normal

MAX_POINTS = 10000

def grow_axon(x, y):
    axon_length = np.random.rayleigh(sigma_l)
    angle = random.uniform(0, 2 * np.pi)

    axon_points = np.zeros((MAX_POINTS, 2))
    axon_points[0] = [x, y]
    count = 1
    total_length = 0

    while count < MAX_POINTS and total_length < axon_length:
        angle += np.random.normal(0, sigma_theta)
        last_point = axon_points[count - 1]
        next_step = segment_length * np.array([np.cos(angle), np.sin(angle)])
        new_point = last_point + next_step

        if np.linalg.norm(new_point) > circle_radius - margin:
            normal = new_point / np.linalg.norm(new_point)
            reflected_direction = reflect_direction(next_step, normal)
            new_point = last_point + segment_length * reflected_direction
            angle = np.arctan2(reflected_direction[1], reflected_direction[0])

        actual_step = np.linalg.norm(new_point - last_point)
        total_length += actual_step

        axon_points[count] = new_point
        count += 1

    return axon_points[:count], total_length

def does_axon_intersect_dendrite(axon_path, neuron_center, dendritic_radius):
    for (x, y) in axon_path:
        if np.linalg.norm(np.array([x, y]) - np.array(neuron_center)) <= dendritic_radius:
            return True
    return False

def analyze_network(A):
    G = nx.DiGraph(A)
    density = nx.density(G)
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())
    avg_in_degree = np.mean(list(in_degrees.values()))
    avg_out_degree = np.mean(list(out_degrees.values()))

    avg_path_length = None
    if nx.is_strongly_connected(G):
        avg_path_length = nx.average_shortest_path_length(G)

    clustering_coeff = nx.average_clustering(G)
    return density, avg_in_degree, avg_out_degree, avg_path_length, clustering_coeff

def create_adjacency_matrix(neuron_connections, N):
    A = np.zeros((N, N))
    for i, connected_neurons in neuron_connections.items():
        for j in connected_neurons:
            A[j, i] = 1   # column = presynaptic, row = postsynaptic
    return A

def find_non_connected_neurons(A):
    N = A.shape[0]
    non_connected = []
    for i in range(N):
        # column i = outgoing from i, row i = incoming to i
        has_outgoing = np.any(A[:, i])
        has_incoming = np.any(A[i, :])
        if not has_outgoing and not has_incoming:
            non_connected.append(i)
    return non_connected

# ---------- Step 1: neuron positions ----------
neuron_positions = []
while len(neuron_positions) < num_neurons:
    x = random.uniform(-circle_radius + neuron_radius, circle_radius - neuron_radius)
    y = random.uniform(-circle_radius + neuron_radius, circle_radius - neuron_radius)
    if is_valid_position(x, y, neuron_positions):
        neuron_positions.append((x, y))
neuron_positions = sorted(neuron_positions, key=lambda pos: pos[0])


with open("position_N_%.5d_alpha_%.2f_seed_%.4d_diam_%.1lf.dat"
          % (num_neurons, alpha, seed, circle_diameter), "w") as f:
    f.write("# Neuron_index\tx_position_mm\ty_position_mm\n")
    for i, (x, y) in enumerate(neuron_positions):
        f.write(f"{i}\t{x:.6f}\t{y:.6f}\n")

# ---------- Step 2: dendritic radii + axon growth ----------
dendritic_radii = []
axon_lengths = []
axon_paths = []

for (x, y) in neuron_positions:
    dendritic_diameter = np.random.normal(mu, sigma)
    dendritic_radii.append(dendritic_diameter / 2)
    axon_path, length = grow_axon(x, y)
    axon_lengths.append(length)
    axon_paths.append(axon_path)

avg_axon_length = np.mean(axon_lengths)
median_axon_length = np.median(axon_lengths)
std_axon_length = np.std(axon_lengths)

# ---------- Step 3: connections (all-to-all, probability alpha) ----------
connections = {i: [] for i in range(num_neurons)}

for i in range(num_neurons):
    axon = axon_paths[i]
    for j in range(num_neurons):
        if i == j:
            continue
        dendrite_center = np.array(neuron_positions[j])
        dendrite_radius = dendritic_radii[j]

        # Only if the axon physically reaches the dendritic tree
        if does_axon_intersect_dendrite(axon, dendrite_center, dendrite_radius):
            # ... then form the synapse with probability alpha
            if np.random.rand() < alpha:
                connections[i].append(j)

N = num_neurons
A = create_adjacency_matrix(connections, N)
density, avg_in, avg_out, avg_path, clust = analyze_network(A)

# ---------- Step 4: statistics ----------
print("\n--- Network Statistics ---\n")
print(f"Neuronal Density: {neuron_density}")
print(f"Density: {density}")
print(f"Average In-Degree: {avg_in}")
print(f"Average Out-Degree: {avg_out}")
print(f"Average Path Length: {avg_path if avg_path is not None else 'Not all neurons are connected'}")
print(f"Clustering Coefficient: {clust}")
print("Number of Connections: %d" % (np.sum(A)))
print("Average Connectivity: %lf" % (np.sum(A) / num_neurons))

print("\n--- Axon Length Statistics ---\n")
print(f"Mean axon length: {avg_axon_length:.3f} mm")
print(f"Median axon length: {median_axon_length:.3f} mm")
print(f"Standard deviation: {std_axon_length:.3f} mm")
print(f"Target mean (Rayleigh): {mean_length:.3f} mm")
print(f"Target sigma (Rayleigh): {sigma_l:.3f} mm")

with open("network_analysis_N_%.5d_alpha_%.2f_seed_%.4d_diam_%.1lf.dat"
          % (N, alpha, seed, circle_diameter), "w") as f:
    f.write(f"Neuronal Density: {neuron_density}\n")
    f.write(f"Neurons: {num_neurons}\n")
    f.write(f"Culture Area: {culture_area}\n")
    f.write(f"Connection Probability alpha: {alpha}\n")
    f.write(f"Connection Density: {density}\n")
    f.write(f"Average In-Degree: {avg_in}\n")
    f.write(f"Average Out-Degree: {avg_out}\n")
    f.write(f"Average Path Length: {avg_path if avg_path is not None else 'Not all neurons are connected'}\n")
    f.write(f"Clustering Coefficient: {clust}\n")
    f.write("Number of Connections: %d\n" % (np.sum(A)))
    f.write("\n--- Axon Length Statistics ---\n")
    f.write(f"Mean axon length: {avg_axon_length:.3f} mm\n")
    f.write(f"Median axon length: {median_axon_length:.3f} mm\n")
    f.write(f"Standard deviation: {std_axon_length:.3f} mm\n")
    f.write(f"Target mean (Rayleigh): {mean_length:.3f} mm\n")
    f.write(f"Target sigma (Rayleigh): {sigma_l:.3f} mm\n")

# ---------- Step 5: adjacency matrix file ----------
vector = []
for i in range(N):
    for j in range(N):
        if A[i, j] == 1:
            vector.append(i * N + j)

np.savetxt("culture_matrix_N_%.5d_alpha_%.2f_seed_%.4d_diam_%.1lf.dat"
           % (N, alpha, seed, circle_diameter),
           vector, delimiter=' ', newline="\n", fmt='%d')

isolated_neurons = find_non_connected_neurons(A)

if isolated_neurons:
    print(f"\nNon-connected neurons: {isolated_neurons}")
else:
    print("\nAll neurons are connected to at least one other neuron.")

# ---------- Step 6: plotting ----------
fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])
ax_0 = fig.add_subplot(gs[0])
ax_1 = fig.add_subplot(gs[1])

# --- Left: culture ---
ax_0.set_title(r"N=%.5d $\rho$ = %.2f neu/mm$^2$ $\bar n=%.3f$ $\alpha=%.2f$ "
               % (N, neuron_density, avg_in, alpha))

outer_circle = plt.Circle((0, 0), circle_radius, color='darkgray',
                          fill=False, linestyle='-', lw=2)
ax_0.add_artist(outer_circle)

# Axons behind
for (x, y), axon_path in zip(neuron_positions, axon_paths):
    if len(axon_path) > 1:
        ax_0.plot(axon_path[:, 0], axon_path[:, 1],
                  color='darkgreen', linewidth=0.5, alpha=0.6, zorder=1)

# Dendritic trees
for (x, y), dendritic_radius in zip(neuron_positions, dendritic_radii):
    ax_0.add_artist(plt.Circle((x, y), dendritic_radius,
                               color='blue', fill=True, alpha=0.1, zorder=2))

# Somata
for (x, y) in neuron_positions:
    ax_0.add_artist(plt.Circle((x, y), neuron_radius,
                               color='darkgreen', fill=True, zorder=3))

ax_0.set_xlim(-circle_radius - .2, circle_radius + .2)
ax_0.set_ylim(-circle_radius - .2, circle_radius + .2)
ax_0.set_aspect('equal', 'box')
ax_0.set_axis_off()

# Scale bar (1 mm)
sb_len = 1.0
sb_x0 = -circle_radius
sb_y0 = -circle_radius - 0.1
ax_0.plot([sb_x0, sb_x0 + sb_len], [sb_y0, sb_y0], 'k-', lw=2)
ax_0.text(sb_x0 + sb_len / 2, sb_y0 - 0.08, f'{sb_len:.1f} mm', ha='center')

# --- Right: adjacency matrix ---
ax_1.set_title("Adjacency Matrix")

A_color = np.ones((N, N, 3))   # white background
for i in range(N):
    for j in range(N):
        if A[i, j] > 0:
            A_color[i, j] = [0, 0.5, 0]   # dark green = connection
        elif i == j and i in isolated_neurons:
            A_color[i, :] = [0.7, 0.7, 0.7]
            A_color[:, j] = [0.7, 0.7, 0.7]

ax_1.imshow(A_color, interpolation="nearest")
ax_1.set_xlabel("Presynaptic Neuron (Axon)")
ax_1.set_ylabel("Postsynaptic Neuron (Dendrite)")

ax_1.set_xticks(np.arange(N))
ax_1.set_yticks(np.arange(N))
ax_1.set_xticklabels(np.arange(N))
ax_1.set_yticklabels(np.arange(N))
plt.setp(ax_1.get_xticklabels(), rotation=90)

legend_elements_adj = [
    Patch(facecolor='darkgreen', label='Connection'),
    Patch(facecolor='gray', alpha=0.5, label='Isolated neurons'),
]
ax_1.legend(handles=legend_elements_adj, loc='upper right', fontsize=8)

plt.subplots_adjust(left=None, bottom=.1, right=.9, top=None, wspace=.25, hspace=.25)
width, height = 28, 24
fig.set_size_inches(width / 2.54, height / 2.54)
plt.savefig('fig_culture_N_%.5d_alpha_%.2f_seed_%.4d.png'
            % (N, alpha, seed), dpi=256, bbox_inches='tight')
plt.close(fig)