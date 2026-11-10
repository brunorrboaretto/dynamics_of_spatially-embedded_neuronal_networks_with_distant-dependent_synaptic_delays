#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#define PI acos(-1)

#define N 500
#define Neq 5
#define NN N * Neq

double gaussrand();

void eqs(double t,double *Ydot, double *Y,double h, int *A, int num_connections, double eps_A, double **Yr);

int main(int argc, char *argv[])
{
  int seed1 = atoi (argv[1]);
  int seed2 = atoi (argv[2]);
  
  int i,j,k;
  double tempo =  10000;
  double trans =  2000;
  double h=0.01;
  int    H=(int)(1/h);
  int Tempo=tempo*H;
  int Trans=trans*H;
  
  int t;
  char *filename_1 = malloc(sizeof(char *) * (200));
  char *filename_X = malloc(sizeof(char *) * (200));

  //Free parameters

  //double eps_A= strtod(argv[3], NULL);
  double D =  1.5;//strtod(argv[3], NULL);
  double vel = strtod(argv[3],NULL);
  double max_distance = strtod(argv[4],NULL);

  double **Yr1 = (double **)malloc(N * sizeof(double *));
  double **Yr2 = (double **)malloc(N * sizeof(double *));
  for (i = 0; i < N; i++)
    {
      Yr1[i] = (double *)malloc(N * sizeof(double));
      Yr2[i] = (double *)malloc(N * sizeof(double));
    }
  int **T_delay = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    T_delay[i] = (int *)malloc(N * sizeof(int));
  
  double eps_A;
  double delta_eps=0.025;
  int E;
  printf("AAA\n");

  sprintf(filename_1, "HH_ddd_cplx_cont_SYNC_N_%.5d_D_%.3lf_vel_%s_diam_%.2lf_seed_%.4d.dat",N,D,argv[3],max_distance,seed1);
  
  FILE *p5 = fopen(filename_1, "w");

  char filename[500];

  int value;
  int num_connections = 0;
  int *A;

  snprintf(filename,sizeof(filename),
	   "culture_matrix_N_%05d_alpha_0.10_seed_%04d_diam_%.1lf.dat",N,seed2,max_distance);
  FILE *fp = fopen(filename, "r");
  printf("%s\n", filename);

  if(fp == NULL)
    {
      printf("Erro ao abrir o arquivo.\n");
      return 1;
    }
  
  while(fscanf(fp,"%d",&value) == 1)
    num_connections++;
  
  printf("Numero de elementos = %d\n",num_connections);
  
  A = (int *) malloc(num_connections*sizeof(int));
  rewind(fp);
  for(i=0;i<num_connections;i++)
    fscanf(fp,"%d",&A[i]);

  fclose(fp);

  char position_filename[500];
 
  double *x_pos = (double *)malloc(N * sizeof(double));
  double *y_pos = (double *)malloc(N * sizeof(double));

  snprintf(position_filename, sizeof(position_filename),
	   "position_N_%05d_alpha_0.10_seed_%04d_diam_%.1lf.dat",
	   N, seed2,max_distance);
  
  FILE *fp_pos = fopen(position_filename, "r");
  
  if (fp_pos == NULL)
    {
      printf("Erro ao abrir arquivo de posições: %s\n", position_filename);
      return 1;
    }

  printf("%s\n", position_filename);
  /* Pula a primeira linha:
     # Neuron_index x_position_mm y_position_mm
  */
  char line[300];
  fgets(line, sizeof(line), fp_pos);

  int neuron_index;

  for (i = 0; i < N; i++)
    {
      if (fscanf(fp_pos, "%d %lf %lf",
		 &neuron_index, &x_pos[i], &y_pos[i]) != 3)
	{
	  printf("Erro ao ler posição do neurônio %d\n", i);
	  fclose(fp_pos);
	  return 1;
	}
    }

  fclose(fp_pos);

  for (k = 0; k < num_connections; k++)
    {
      i = A[k] / N;
      j = A[k] % N;

      double dx = x_pos[i] - x_pos[j];
      double dy = y_pos[i] - y_pos[j];

      double distance = sqrt(dx * dx + dy * dy);

      double tau_ij = distance / vel;

      T_delay[i][j] = (int)(tau_ij / h);
    }


  
  double **r_d = (double **)malloc(N * sizeof(double *));
  double **r_da = (double **)malloc(N * sizeof(double *));
  for (i = 0; i < N; i++)
    {
      r_d[i] = (double *)malloc(Tempo * sizeof(double));
      r_da[i] = (double *)malloc(Tempo * sizeof(double));
    }
  
  double *Ydot =(double*)malloc(NN*sizeof(double)); 
  double *Yaux =(double*)malloc(NN*sizeof(double)); 
  double *Y    =(double*)malloc(NN*sizeof(double));
  double *Yc   =(double*)malloc(NN*sizeof(double));
  double *K1   =(double*)malloc(NN*sizeof(double)); 
  double *K2   =(double*)malloc(NN*sizeof(double));
  double *K3   =(double*)malloc(NN*sizeof(double)); 
  double *K4   =(double*)malloc(NN*sizeof(double)); 

  double *L1   =(double*)malloc(NN*sizeof(double)); 
  double *X1 = (double *)malloc(N * sizeof(double));
  double *X2 = (double *)malloc(N * sizeof(double));
 
  int T_delay_min_int = INT_MAX;
  int T_delay_max_int = 0;
  
  
  for (k = 0; k < num_connections; k++)
    {
      i = A[k] / N;
      j = A[k] % N;

      if (T_delay[i][j] < T_delay_min_int)
	T_delay_min_int = T_delay[i][j];

      if (T_delay[i][j] > T_delay_max_int)
	T_delay_max_int = T_delay[i][j];
     
    }

  printf("Minimum integer delay = %d steps\n", T_delay_min_int);
  printf("Maximum integer delay = %d steps\n", T_delay_max_int);

  //ICs----------------------------------------//
  
  srand(seed1);
  
  for (i = 0; i < NN; i += Neq)
    {
      Yc[i+0] = -80.0 + ((double)(rand() % 1000000) / 1000000.0) * 80.0;
      Yc[i+1] = 0.0 + ((double)(rand() % 1000000) / 1000000.0);
      Yc[i+2] = 0.0 + ((double)(rand() % 1000000) / 1000000.0);
      Yc[i+3] = 0.0 + ((double)(rand() % 1000000) / 1000000.0);
      Yc[i+4] = 0.0 + ((double)(rand() % 1000000) / 1000000.0);
    }
  
  int nspikes = 10000;
  int **k1 = (int **)malloc(N*sizeof(int *));
  for(i=0;i<N;i++)
    {
      k1[i]=(int *)malloc(nspikes*sizeof(int));
    }
  int l[N];
  double *COS = (double *)malloc(Tempo * sizeof(double));
  double *SIN = (double *)malloc(Tempo * sizeof(double));
  double *COS2 = (double *)malloc(Tempo * sizeof(double));
  double *SIN2 = (double *)malloc(Tempo * sizeof(double));
  double *COS3 = (double *)malloc(Tempo * sizeof(double));
  double *SIN3 = (double *)malloc(Tempo * sizeof(double));
  double *COS4 = (double *)malloc(Tempo * sizeof(double));
  double *SIN4 = (double *)malloc(Tempo * sizeof(double));

  double *fase = (double *)malloc(Tempo * sizeof(double));
  
  //double **fase = malloc(N * sizeof(double *));
  //for (i = 0; i < N; i++)
  //  fase[i] = malloc(Tempo * sizeof(double));

  double spike_CV = 0.0;
  double xi = 0;
  double cmed = 0.0;
  double max_cmed = -10000.0;
  double min_cmed = 10000.0;

  for (i = 0; i < N; i++)
    {
      l[i] = 0;
      for(j=0;j<N;j++)
	{
	  Yr1[i][j] = 0.;
	  Yr2[i][j] = 0.;
	}
      
      for (t = 0; t < Tempo; t++)
	{
	  r_d[i][t] = 0.;
	  r_da[i][t]= 0.;
	}
      for(t = 0;t<nspikes;t++)
	{
	  k1[i][t]=0;
	}
    }

  for (i = 0; i < NN; i++)
    {
      Y[i] = Yc[i];
      L1[i]=0.0;
    }
  
  for (t = 0; t < Tempo; t++) //primeira iterada
    {
      COS[t] = 0.;
      SIN[t] = 0.;
      COS2[t] = 0.;
      SIN2[t]=  0.;
      COS3[t] = 0.;
      SIN3[t]=  0.;
      COS4[t] = 0.;
      SIN4[t]=  0.;
      fase[t]=0.0;
      for(i=0;i<N;i++)
	{
	  xi = gaussrand();
	  L1[i*Neq]=sqrt(h)*xi*D;
	  //fase[i][t]=0.0;
	}
      
      eqs(t * h + 0.0 * h, K1, Y   , h,A,num_connections,0.0,Yr1);
      for (i = 0; i < NN; i++) Yaux[i] = Y[i] + h*K1[i]+L1[i]; //aqui entra o l
      eqs(t * h + 1.0 * h, K2, Yaux, h,A,num_connections,0.0,Yr2);

      for (i = 0; i < NN; i++)
	{
	  Ydot[i] = Y[i] + (h / 2.0) * (K1[i] + K2[i])+L1[i];
	  if(i%Neq==0)
	    {
	      X2[i / Neq] = Ydot[i];
	      X1[i / Neq] = Y[i];
	    }
	  Y[i] = Ydot[i];
	}
      
      for (i = 0; i < N; i++)
	{
	  r_da[i][t] = Y[i * Neq + 4];
	}
    }
  
  for(E=0;E<=40;E++)
    {
      eps_A = E*delta_eps; 
      
      max_cmed = -10000.0;
      min_cmed = 10000.0;
      //sprintf(filename_X, "new_results_NC/HH_ddd_cplx_cont_CORR_N_%.5d_D_%.3lf_vel_%s_diam_%.2lf_seed_%.4d_eps_%.3lf.dat",N,D,argv[3],max_distance,seed1,eps_A);
      //FILE *pX = fopen(filename_X, "w");
      
      for(t=0;t<Tempo;t++)
	{
	  cmed = 0.0;
	  
	  for (k = 0; k < num_connections; k++)
	    {
	      i = A[k] / N;
	      j = A[k] % N;

	      int td = T_delay[i][j];

	      if (td == 0)
		{
		  Yr1[i][j] = Y[j * Neq + 4];
		}
	      else if (t >= td)
		{
		  Yr1[i][j] = r_d[j][t - td];
		}
	      else
		{
		  Yr1[i][j] = r_da[j][Tempo - td + t];
		}
	    }
	  
	  for(i=0;i<N;i++)
	    {
	      xi = gaussrand();
	      L1[i*Neq]=sqrt(h)*xi*D;
	    }
      
	  eqs(t * h + 0.0 * h, K1, Y   , h,A,num_connections, eps_A,Yr1);
	  for (i = 0; i < NN; i++) Yaux[i] = Y[i] + h*K1[i]+L1[i]; //aqui entra o l
	  
	  for (k = 0; k < num_connections; k++)
	    {
	      i = A[k] / N;
	      j = A[k] % N;
	      
	      int td = T_delay[i][j];
	      
	      if (td == 0)
		{
		  Yr2[i][j] = Yaux[j * Neq + 4];
		}
	      else if (t + 1 >= td)
		{
		  Yr2[i][j] = r_d[j][t + 1 - td];
		}
	      else
		{
		  Yr2[i][j] = r_da[j][Tempo - td + t + 1];
		}
	    }
	  eqs(t * h + 1.0 * h, K2, Yaux, h,A,num_connections,eps_A,Yr2);

	  for (i = 0; i < NN; i++)
	    {
	      Ydot[i] = Y[i] + (h / 2.0) * (K1[i] + K2[i])+L1[i];
	      if(i%Neq==0)
		{
		  X2[i / Neq] = Ydot[i];
		  X1[i / Neq] = Y[i];
		}
	      Y[i] = Ydot[i];
	    }

	  for (i = 0; i < N; i++)
	    {
	      if(X2[i]>-20.0 && X1[i]<-20.0)
		{
		  k1[i][l[i]] = t;
		  l[i]++;
		}

	      r_d[i][t] = Y[i * Neq + 4];

	      cmed+=Y[i*Neq];
	    }

	  cmed/=N;

	  if(cmed>max_cmed)
	    max_cmed=cmed;
	  if(cmed<min_cmed)
	    min_cmed=cmed;   
	}
      printf("AAAA\n");

      int tmax = Tempo;
      int tmin = Trans;

      for(i=0;i<N;i++)
	{
	  if(k1[i][l[i]-1] < tmax)
	    tmax = k1[i][l[i]-1];
	  if(k1[i][0] > tmin)
	    tmin = k1[i][0];
	}
      printf("%lf %d %d %d\n",eps_A,tmax,tmin,l[0]);
      
      double soma_ISI = 0.0;
      long int n_ISI = 0;

      int L;
      for(i=0;i<N;i++)
	{
	  for(L=0;L<l[i]-1;L++)
	    {
	      int k_0 = k1[i][L];
	      int k_1 = k1[i][L+1];
	      
	      soma_ISI += (k_1-k_0)*h;
	      n_ISI++;
	      
	      for(k=k_0;k<k_1;k++)
		{
		  fase[k]=2*PI*L + 2*PI*(k-k_0)/(k_1-k_0);
		  
		  COS[k]+=cos(fase[k]);
		  SIN[k]+=sin(fase[k]);
		  COS2[k]+=cos(2*fase[k]);
		  SIN2[k]+=sin(2*fase[k]);
		  COS3[k]+=cos(3*fase[k]);
		  SIN3[k]+=sin(3*fase[k]);
		  COS4[k]+=cos(4*fase[k]);
		  SIN4[k]+=sin(4*fase[k]);
		  
		  //fase[i][k]=2*PI*L + 2*PI*(k-k_0)/(k_1-k_0);
		  
		  //COS[k]+=cos(fase[i][k]);
		  //SIN[k]+=sin(fase[i][k]);
		  //COS2[k]+=cos(2*fase[i][k]);
		  //SIN2[k]+=sin(2*fase[i][k]);
		  //COS3[k]+=cos(3*fase[i][k]);
		  //SIN3[k]+=sin(3*fase[i][k]);
		  //COS4[k]+=cos(4*fase[i][k]);
		  //SIN4[k]+=sin(4*fase[i][k]);
		}
	    }
	}

      double T_medio = soma_ISI/n_ISI;
      double f_media = 1.0/T_medio;
      
      double R;
      double R2;
      double R3;
      double R4;
      
      double rmed=0.0;
      double rmed2=0.0;
      double rmed3=0.0;
      double rmed4=0.0;

      //int n_samples = 0;
      //double **C = malloc(N * sizeof(double *));
      //double **S = malloc(N * sizeof(double *));

      //for (i = 0; i < N; i++)
      //{
      //	  C[i] = calloc(N, sizeof(double));
      //	  S[i] = calloc(N, sizeof(double));
      //	}
      
      for (k = tmin; k < tmax; k++)
	{
	  R = sqrt((pow(COS[k], 2) + pow(SIN[k], 2)))/N;
	  R2 = sqrt((pow(COS2[k], 2) + pow(SIN2[k], 2)))/N;
	  R3 = sqrt((pow(COS3[k], 2) + pow(SIN3[k], 2)))/N;
	  R4 = sqrt((pow(COS4[k], 2) + pow(SIN4[k], 2)))/N;
	  /*
	  if(k%H==0)
	    {
	      n_samples++;
	      for(i=0;i<N;i++)
		for (j = i + 1; j < N; j++)
		  {
		    C[i][j] += cos(fase[i][k] - fase[j][k]);
		    S[i][j] += sin(fase[i][k] - fase[j][k]);
		  }
	    }
	  */
	  rmed+=R;
	  rmed2+=R2;
	  rmed3+=R3;
	  rmed4+=R4;

	  //if(k%H==0)
	  //fprintf(p4,"%lf %lf %lf %lf %lf\n",tempo*E+k*h,R,R2,R3,R4);
	}

      rmed/=(tmax-tmin);
      rmed2/=(tmax-tmin);
      rmed3/=(tmax-tmin);
      rmed4/=(tmax-tmin);

      
      fprintf(p5,"%lf %lf %lf %lf %lf \n",eps_A,fabs(max_cmed-min_cmed),rmed,rmed2,T_medio);	  
      printf("%lf %lf %lf %lf %lf %lf \n",t*h,eps_A,fabs(max_cmed-min_cmed),rmed,rmed2,T_medio);
      fflush(p5);
      /*      
      for(i=0;i<N;i++)
	for (j = i + 1; j < N; j++)
	  {
	    C[i][j] /= (double)n_samples;
	    S[i][j] /= (double)n_samples;
	  }
      
      double delta_d = 0.01;   // largura do bin em mm
      

      int n_bins = (int)(max_distance / delta_d) + 1;

      double *C_dist = calloc(n_bins, sizeof(double));
      double *S_dist = calloc(n_bins, sizeof(double));
      int *N_dist = calloc(n_bins, sizeof(int));

      for (i = 0; i < N; i++)
	{
	  for (j = i + 1; j < N; j++)
	    {
	      double dx = x_pos[i] - x_pos[j];
	      double dy = y_pos[i] - y_pos[j];

	      double distance = sqrt(dx*dx + dy*dy);

	      int bin = (int)(distance / delta_d);
      
	      if (bin < 0 || bin >= n_bins)
		{
		  fprintf(stderr,
			  "Erro: distancia fora dos bins: i=%d j=%d "
			  "distance=%.10f bin=%d n_bins=%d\n",
			  i, j, distance, bin, n_bins);
		  exit(EXIT_FAILURE);
		}
	      C_dist[bin] += C[i][j];
	      S_dist[bin] += S[i][j];
	      N_dist[bin]++;
	    }
	}
      
      for (i = 0; i < n_bins; i++)
	{
	  if (N_dist[i] > 0)
	    {
	      C_dist[i] /= N_dist[i];
	      S_dist[i] /= N_dist[i];
	    }
	}

      for (i = 0; i < n_bins; i++)
	{
	  if (N_dist[i] > 0)
	    {
	      double distance_center = (i + 0.5) * delta_d;
	      
	      fprintf(pX, "%lf %lf %lf %lf %d\n",
		      distance_center,
		      C_dist[i],S_dist[i], sqrt(C_dist[i]*C_dist[i] + S_dist[i]*S_dist[i]),N_dist[i]);
	    }
	}
      */
      //fflush(pX);
      //fclose(pX);

      //for (i = 0; i < N; i++)
      //	free(C[i]);

      //free(C);
      //free(C_dist);
      //free(N_dist);

      
      for(i=0;i<N;i++)
	{
	  l[i] = 0;
	  for(j=0;j<N;j++)
	    {
	      Yr1[i][j] = 0.;
	      Yr2[i][j] = 0.;
	    }
	  for (t = 0; t < Tempo; t++)
	    {
	      r_da[i][t]= r_d[i][t];
	      r_d[i][t] = 0.;
	      COS[t]=0.;
	      SIN[t]=0.;
	      COS2[t]=0.;
	      SIN2[t]=0.;
	      COS3[t]=0.;
	      SIN3[t]=0.;
	      COS4[t]=0.;
	      SIN4[t]=0.;
	    }
	  
	  for(t = 0;t<nspikes;t++)
	    {
	      k1[i][t]=0;
	    }
	}
    }

  /* Fecha os arquivos de saida */
  fclose(p5);

  /* Libera as matrizes dos atrasos */
  for (i = 0; i < N; i++)
    {
      free(Yr1[i]);
      free(Yr2[i]);
      free(T_delay[i]);
    }

  free(Yr1);
  free(Yr2);
  free(T_delay);

  /* Libera o historico da variavel neuronal */
  for (i = 0; i < N; i++)
    {
      free(r_d[i]);
      free(r_da[i]);
    }

  free(r_d);
  free(r_da);

  /* Libera os vetores do integrador */
  free(Ydot);
  free(Yaux);
  free(Y);
  free(Yc);
  free(K1);
  free(K2);
  free(K3);
  free(K4);
  free(L1);

  free(X1);
  free(X2);

  /* Libera os dados das conexoes e posicoes */
  free(A);
  free(x_pos);
  free(y_pos);

  /* Libera os tempos de disparo */
  for (i = 0; i < N; i++)
    free(k1[i]);

  free(k1);

  /* Libera os vetores dos parametros de ordem */
  free(COS);
  free(SIN);
  free(COS2);
  free(SIN2);
  free(COS3);
  free(SIN3);
  free(COS4);
  free(SIN4);
  
  /* Libera a matriz de fases */
  //for (i = 0; i < N; i++)
  //  free(fase[i]);

  free(fase);

  /* Libera os nomes dos arquivos */
  free(filename_1);
  free(filename_X);

  return 0;
  
}

int **create_connections(double p, int seed, int *num_connections) {
  int **A;
  *num_connections = 0;
  srand(seed);

  int i,j;
  
  // Allocate NxN matrix
  A = (int **)malloc(N * sizeof(int *));
  if (A == NULL)
    {
      fprintf(stderr, "Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }

  for (i = 0; i < N; i++)
    {
      A[i] = (int *)malloc(N * sizeof(int));
      if (A[i] == NULL)
	{
	  fprintf(stderr, "Memory allocation failed.\n");
	  exit(EXIT_FAILURE);
	}
    }
  
  // Fill matrix with random directed connections
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
	{
	  if (i == j)
	    A[i][j] = 0; // avoid self-connection
	  else
	    {
	      double r = (double)rand() / RAND_MAX;
	      if (r < p)
		{
		  A[i][j] = 1;
		  (*num_connections)++;
		}
	      else
		A[i][j] = 0;
	    }
	}
    }
  
  // Check for isolated neurons
  for (i = 0; i < N; i++)
    {
      int in_deg = 0, out_deg = 0;
      for (j = 0; j < N; j++)
	{
	  out_deg += A[i][j];
	  in_deg += A[j][i];
	}
      if (in_deg == 0 && out_deg == 0)
	{
	  printf("Warning: neuron %d is isolated.\n", i);
	}
    }
  
  return A;
}

int *matrix_to_vector(int **A, int *vector_size)
{
  int count = 0;
  int i,j;
  // First, count how many connections (1s) exist
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (A[i][j] == 1)
	count++;
	
  // Allocate memory for the vector
  int *connections = (int *)malloc(count * sizeof(int));
  if (connections == NULL)
    {
      fprintf(stderr, "Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }
  
  // Fill vector with positions of 1s
  int idx = 0;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (A[i][j] == 1) 
	connections[idx++] = i * N + j;
 
  *vector_size = count;
  return connections;
}



double gaussrand()
{
  static double U, V;
  static int phase = 0;
  double Z;

  if(phase == 0) {
    U = (rand() + 1.) / (RAND_MAX + 2.);
    V = rand() / (RAND_MAX + 1.);
    Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
  } else
    Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

  phase = 1 - phase;

  return Z;
}

void eqs(double t,double *Ydot, double *Y,double h,  int *A, int num_connections, double eps_A,double **Yr)
{
  
  int i,j,k;
  int count;
  
  //potenciais de Nernst de cada íon em mV
  double E_K  =-77.0;
  double E_Na = 50.0;
  double E_vaz=-54.4;
  
  //codutancias de cada íon
  double g_K  = 36.0;
  double g_Na =120.0;
  double g_vaz=  0.3;
  
  //correntes de cada íon (associados a uma edo)
  
  double a_n,a_m,a_h;
  double b_n,b_m,b_h;
  double e_Vsyn = 40.0;
  double i_Vsyn = -80.0;
  
  double Tr = 0.4;
  double Td = 2.0;
  
  double I = 0.0;//10.0;

  double IA[N];

  for (i = 0; i < N; i++)
    IA[i] = 0.0;
  
  for(k=0;k<num_connections;k++)
    {
      i = A[k]/N;
      j = A[k]%N;
      
      IA[i] +=Yr[i][j];
    }
  
  double Norm = (float)(num_connections/N);

  for(i=0;i<NN;i+=Neq)
    {  
      a_n=0.01*(Y[i+0]+55.0)/(1-exp(-(Y[i+0]+55)/10.0));
      a_m=0.10*(Y[i+0]+40.0)/(1-exp(-(Y[i+0]+40)/10.0));
      a_h=0.07*exp(-(Y[i+0]+65.0)/20.0);
      
      b_n=0.125*exp(-(Y[i+0]+65.0)/80.0);
      b_m=4.000*exp(-(Y[i+0]+65.0)/18.0);
      b_h=1/(1+exp(-(Y[i+0]+35.0)/10.0));
       
      Ydot[i+0] =-g_Na*pow(Y[i+2],3)*Y[i+3]*(Y[i+0]-E_Na)-g_K*pow(Y[i+1],4)*(Y[i+0]-E_K)-g_vaz*(Y[i+0]-E_vaz) + (eps_A/(Norm))*(e_Vsyn-Y[i])*(IA[i/Neq]);
      Ydot[i+1] = a_n*(1-Y[i+1])-b_n*Y[i+1];
      Ydot[i+2] = a_m*(1-Y[i+2])-b_m*Y[i+2];
      Ydot[i+3] = a_h*(1-Y[i+3])-b_h*Y[i+3];
      Ydot[i+4] = ((1 / Tr) - (1 / Td)) * (1 - Y[i + 4]) / (1 + exp(-(Y[i] + 20))) - Y[i + 4] / Td;
    }
}


int **create_WS_connections(int K, double p_rewire, int seed, 
			    int *num_connections_total, int *num_regular_connections, int *num_random_connections) {
  int **A;
  *num_connections_total = 0;
  *num_regular_connections = 0;
  *num_random_connections = 0;
  srand(seed);

  int i, j;
  int offset, c, try_node;
    
  if (K >= N) {
    fprintf(stderr, "K must be smaller than N.\n");
    exit(EXIT_FAILURE);
  }

  if (K % 2 != 0) {
    fprintf(stderr, "Warning: K should be even for symmetric connections. Using K=%d\n", K);
  }

  // Allocate NxN matrix
  A = (int **)malloc(N * sizeof(int *));
  if (A == NULL) {
    fprintf(stderr, "Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < N; i++) {
    A[i] = (int *)calloc(N, sizeof(int));
    if (A[i] == NULL) {
      fprintf(stderr, "Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }
  }
    
  // Track whether a connection is regular (0) or random (1)
  int **connection_type = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++) {
    connection_type[i] = (int *)calloc(N, sizeof(int));
  }
    
  // Step 1: Create regular ring lattice
  // Connect each node to its K nearest neighbors (K/2 on each side)
  for (i = 0; i < N; i++) {
    // Connect to K/2 neighbors on the right side
    for (offset = 1; offset <= K/2; offset++) {
      j = (i + offset) % N;
      A[i][j] = 1;
      connection_type[i][j] = 0;  // 0 = regular connection
    }
    // Connect to K/2 neighbors on the left side
    for (offset = 1; offset <= K/2; offset++) {
      j = (i - offset + N) % N;
      A[i][j] = 1;
      connection_type[i][j] = 0;  // 0 = regular connection
    }
  }
    
  // Count initial regular connections
  *num_regular_connections = N * K;
  *num_connections_total = *num_regular_connections;
    
  // Step 2: Rewire edges with probability p_rewire
  for (i = 0; i < N; i++) {
    // First, collect all outgoing connections from node i
    int *outgoing_connections = malloc(K * sizeof(int));
    int conn_count = 0;
        
    // Find all K outgoing connections from this node
    for (j = 0; j < N; j++) {
      if (A[i][j] == 1) {
	outgoing_connections[conn_count++] = j;
      }
    }
        
    // Now rewire each of these connections
    for (c = 0; c < conn_count; c++) {
      int j_old = outgoing_connections[c];
            
      double r = (double)rand() / RAND_MAX;
      if (r < p_rewire) {
	// This connection will become random
	connection_type[i][j_old] = 1;  // Mark as random (even though we'll remove it)
                
	// Remove old connection
	A[i][j_old] = 0;
	(*num_regular_connections)--;
                
	// Find new random target (can be ANY node except self, including old neighbors!)
	int j_new;
	int attempts = 0;
	int max_attempts = 100;
                
	do {
	  j_new = rand() % N;
	  attempts++;
	} while (j_new == i && attempts < max_attempts);
                
	// Special case: p_rewire = 1.0, we might try to rewire to same node
	if (j_new == i) {
	  j_new = (j_new + 1) % N;  // Simple fallback
	}
                
	// Check if connection already exists
	// If it does, we need to find another target
	if (A[i][j_new] == 1) {
	  // Connection already exists, find a different one
	  for (try_node = 0; try_node < N; try_node++) {
	    if (try_node != i && A[i][try_node] == 0) {
	      j_new = try_node;
	      break;
	    }
	  }
	}
                
	// Add new connection as random type
	A[i][j_new] = 1;
	connection_type[i][j_new] = 1;  // 1 = random connection
	(*num_random_connections)++;
      }
    }
        
    free(outgoing_connections);
  }
    
  // Recalculate total connections and verify counts
  *num_connections_total = 0;
  int count_regular = 0;
  int count_random = 0;
    
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (A[i][j] == 1) {
	(*num_connections_total)++;
	if (connection_type[i][j] == 0) {
	  count_regular++;
	} else {
	  count_random++;
	}
      }
    }
  }
    
  // Verify counts match
  if (count_regular != *num_regular_connections || 
      count_random != *num_random_connections ||
      *num_connections_total != *num_regular_connections + *num_random_connections) {
    printf("Warning: Connection counts mismatch!\n");
    printf("Regular: counted=%d, tracked=%d\n", count_regular, *num_regular_connections);
    printf("Random: counted=%d, tracked=%d\n", count_random, *num_random_connections);
    // Update with actual counts
    *num_regular_connections = count_regular;
    *num_random_connections = count_random;
    *num_connections_total = count_regular + count_random;
  }
    
  // Step 3: Check for isolated neurons
  for (i = 0; i < N; i++) {
    int in_deg = 0, out_deg = 0;
    for (j = 0; j < N; j++) {
      out_deg += A[i][j];
      in_deg += A[j][i];
    }
    if (in_deg == 0 && out_deg == 0) {
      printf("Warning: neuron %d is isolated.\n", i);
    }
  }
    
  // Calculate statistics
  double regular_fraction = (double)(*num_regular_connections) / (*num_connections_total);
  double random_fraction = (double)(*num_random_connections) / (*num_connections_total);
    
  printf("\n--- Network Statistics ---\n");
  printf("Total connections: %d\n", *num_connections_total);
  printf("Regular connections: %d (%.2f%%)\n", *num_regular_connections, regular_fraction * 100);
  printf("Random connections: %d (%.2f%%)\n", *num_random_connections, random_fraction * 100);
  printf("Requested rewiring probability: %.2f\n", p_rewire);
  printf("Actual rewired fraction: %.2f\n\n", random_fraction);
    
  // Free connection_type matrix
  for (i = 0; i < N; i++) {
    free(connection_type[i]);
  }
  free(connection_type);

  return A;
}

