#include "matrix.h"

// inner product of 2 vectors
double inner(double* a, double* b, int n);

// normalize vector 
double normalize(double* a, int n);

double power_iteration(Matrix *A, double *v);

double compute_freq(double *params, double ** u, int * number_of_nodes);