#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "lu.h"
#include "eigen.h"
#include "design.h"
#include <gmshc.h>
#include "elasticity.h"
#include <time.h>

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void save_matrix(char *filename, Matrix *A)
{
  FILE *fpt;

  fpt = fopen(filename, "w+");

  for (int i = 0; i < A->m; i++)
  {
    for (int j = 0; j < A->n; j++)
    {
      fprintf(fpt, "%.3e", A->a[i][j]);
      if (j < A->n - 1)
      {
        fprintf(fpt, " ");
      }
    }
    fprintf(fpt, "\n");
  }

  fclose(fpt);
}

Matrix * delete_nodes_from_mat(Matrix *mat, size_t *nodes_list, size_t n_nodes)
{
  int m_new = mat->m - 2 * n_nodes, n_new = mat->n - 2 * n_nodes;
  Matrix *new_mat = allocate_matrix(m_new, n_new);

  int row_offset = 0;
  int i = 0, i_runner = 0;
  while (i < mat->m)
  {
    if (i == 2 * (int)*(nodes_list + row_offset))
    {
      i = i + 2;
      row_offset = (row_offset + 1) % n_nodes;
    }
    else
    {
      int col_offset = 0;
      int j = 0, j_runner = 0;
      while (j < mat->n)
      {
        if (j == 2 * (int)*(nodes_list + col_offset))
        {
          j = j + 2;
          col_offset = (col_offset + 1) % n_nodes;
        }
        else
        {
          new_mat->a[i_runner][j_runner] = mat->a[i][j];
          j++;
          j_runner++;
        }
      }
      i++;
      i_runner++;
    }
  }
  return new_mat;
}

int size_tComparator ( const void * first, const void * second ) {
    size_t firstInt = * (const size_t *) first;
    size_t secondInt = * (const size_t *) second;
    return firstInt - secondInt;
}

// Matrix times vector product
int matrix_times_vector(Matrix *A, double* b, double* res){
  int n = A->n;
  for (int i=0; i<n; i++){
    double b_i = 0.;
    for (int j=0; j<n; j++){
      b_i += A->a[i][j]*b[j];
    }
    res[i] = b_i;
  }
  return 0;
}

// inner product of 2 vectors
double inner(double* a, double* b, int n){
  double res = 0.;
  for (int i=0; i<n; i++) res += a[i]*b[i];
  return res;
}

// normalize vector 
double normalize(double* a, int n){
  double norm = 0.;
  for (int i=0; i<n; i++) norm += a[i]*a[i];
  norm = sqrt(norm);
  for (int i=0; i<n; i++) a[i] /= norm;
  return norm;
}

double power_iteration(Matrix *A, double *v) {
  int n = A->n;

  // initial guess
  double Av[n];
  for (int i=0; i<n; i++) v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v,n);

  double lambda = 1, lambda_prev, diff;
  double rtol = 1e-9; // relative tolerance on lambda

  for(int it = 0; it < 1e4; it++) {
    lambda_prev = lambda;
    // printf("\n===== Iteration %d =====\n", it+1);
    matrix_times_vector(A, v, Av);
    lambda = inner(v, Av, n);
    // printf("λ = %.9e\n", lambda);
    diff = fabs((lambda-lambda_prev) / lambda_prev);
    for(int i = 0; i < n; i++) v[i] = Av[i];
    normalize(v, n);
    if(diff < rtol) break;
  }

  return lambda;
}

/**
 * Builds the tuning fork according to params and computes the frequency
 * */
double compute_freq(double *params, double ** u, int * number_of_nodes)
{
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

  double r1 = params[0];
  double r2 = params[1];
  double e = params[2];
  double l = params[3];
  double mesh_size = 0.1;

  printf("computing frequency using parameters r1: %.3e, r2: %.3e, e: %.3e, l: %.3e\n", r1, r2, e, l);

  designHalfTuningFork(r1, r2, e, l, mesh_size, NULL);

  Matrix *K, *M;
  size_t *boundary_nodes;
  size_t n_boundary_nodes;
  double *coord;
  assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);
  int n_nodes = K->n / 2;

	int * perm = malloc(n_nodes * sizeof(int));
	compute_permutation(perm, coord, n_nodes);

  int * invperm = malloc(n_nodes * sizeof(int));
	for(int i = 0; i < n_nodes; i++) 
		invperm[perm[i]] = i;

  Matrix * K_perm = allocate_matrix(2 * n_nodes, 2 * n_nodes);
  Matrix * M_perm = allocate_matrix(2 * n_nodes, 2 * n_nodes);

  int node_i, node_j;
  for(int i = 0; i < 2* n_nodes; i++) {
    node_i = i/2;
    for(int j = 0; j < 2* n_nodes; j++) {
      node_j = j/2;
      K_perm->a[2*invperm[node_i] + i%2][2*invperm[node_j] + j%2] = K->a[i][j];
      M_perm->a[2*invperm[node_i] + i%2][2*invperm[node_j] + j%2] = M->a[i][j];
    }
  }

  size_t * perm_boundary_nodes = (size_t*) malloc(n_boundary_nodes * sizeof(size_t));
  for(int i = 0; i < n_boundary_nodes; i++) {
    perm_boundary_nodes[i] = invperm[boundary_nodes[i]];
  }

  qsort(perm_boundary_nodes, n_boundary_nodes, sizeof(size_t), size_tComparator);

  Matrix * K_prime = delete_nodes_from_mat(K_perm, perm_boundary_nodes, n_boundary_nodes);
  Matrix * M_prime = delete_nodes_from_mat(M_perm, perm_boundary_nodes, n_boundary_nodes);

  // Get the bandwidth
	int band = 0;
	for (int i = 0; i < K_prime->m; i++) {
		for (int j = 0; j < K_prime->n; j++) {
			if (K_prime->a[i][j] != 0 && abs(j-i) > band) {
				band = abs(j-i);
			}
		}
	}

  BandMatrix * K_band = allocate_band_matrix(K_prime->m, band);
  for (int i = 0; i < K_prime->m; i++) {
		int l1 = MAX(0, i-band);
		int l2 = MIN(K_prime->m-1, i+band);
		for (int j = l1; j <= l2; j++) {
			K_band->a[i][j] = K_prime->a[i][j];
		}
	} 

  Matrix *X = allocate_matrix(K_prime->m, K_prime->n); // X = K^-1M

  lu_band(K_band);

  #pragma omp parallel
  {
    #pragma omp for
    for (int i = 0; i < K_band->m; i++)
    {
      double y[K_band->m];
      for (int j = 0; j < K_band->m; j++)
      {
        y[j] = M_prime->a[j][i];
      }
      solve_band(K_band, y);
      for (int j = 0; j < K_prime->m; j++)
      {
        X->a[j][i] = y[j];
      }
    }
  }

  // Power iteration + deflation to find k largest eigenvalues
  Matrix *A = X; // alias for conciseness
  double *v = malloc(A->m * sizeof(double));
  double lambda, freq;

  lambda = power_iteration(A, v);
  double pi = 3.1415;
  freq = 1. / (2 * pi * sqrt(lambda));

  printf("Frequency computed: %.10f\n", freq);

  if(u) {
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0; 
    for(int i = 0; i < K->m/2; i++) {
      if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        i_bnd++;
        continue;
      }
      vall[2*(perm[i])]   = v[2*iv];
      vall[2*(perm[i])+1] = v[2*iv+1];
      iv++;
    }
    *u = vall;
    *number_of_nodes = n_nodes;
  }

  // free
  free_matrix(K_perm);
  free_matrix(M_perm);
  free_matrix(K);
  free_matrix(M);
  free(v);
  free_matrix(X);
  free_band_matrix(K_band);
  free_matrix(K_prime);
  free_matrix(M_prime);
  free(perm_boundary_nodes);
  free(invperm);
  free(perm);
  free(coord);
  free(boundary_nodes);  

  return freq;
}
