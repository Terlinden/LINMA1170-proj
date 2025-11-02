#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "matrix.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix * mat = (Matrix*) malloc(sizeof(Matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) calloc(m*n, sizeof(double));
	mat->a = (double**) malloc(m*sizeof(double*));
	for(int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
	return mat;
}

BandMatrix * allocate_band_matrix(int m, int k) {
	// To make things simple, we allocate 2k+1 for each row
	// (even though we need less than that for the k first and last rows).
	// We want to have a[i][i] == data[k + i*(2*k+1)]

	BandMatrix * mat = (BandMatrix*) malloc(sizeof(BandMatrix));
	mat->m = m, mat->k = k;
	mat->data = (double*) calloc(m*(2*k+1), sizeof(double));
	mat->a = (double**) malloc(m*sizeof(double*));
	// This is the tricky part :-)
	// We want a[i][i] == data[(2*k+1)*i + k]
	// which means a[i] + i == data + (2*k+1)*i + k
	// and therefore a[i] == data + k + 2*k*i
	for(int i = 0; i < m; i++)
		mat->a[i] = mat->data + k + 2*k*i;
	
	return mat;
}

void free_matrix(Matrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}

void free_band_matrix(BandMatrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}

void print_vector(double * v, int n) {
	for(int i = 0; i < n; i++)
		printf("%.3e ", v[i]);
	printf("\n");
}

void print_matrix(Matrix * A) {
	for(int i = 0; i < A->m; i++)
		print_vector(A->a[i], A->n);
}

int is_symmetric(Matrix * K) {
	int symmetric = 1;
	for(int i = 0; i < K->m; i++)
		for(int j = i+1; j < K->n; j++)
			if(fabs((K->a[i][j] - K->a[j][i]) / K->a[i][j]) > 1e-12) {
				printf("%d %d\n", i, j);
				printf("%lf %lf\n", K->a[i][j], K->a[j][i]);
				symmetric = 0;
			}
	return symmetric;
}

typedef struct {
	int i;		// index
	double x,y;	// coordinates
} Node;

// Comparateur
int cmp(const void * a, const void * b) {
	
	Node * na = (Node *) a;
	Node * nb = (Node *) b;
	if (na->y > nb->y) return 1;
	else return -1;
}

int compute_permutation(int * perm, double * coord, int n_nodes) {

	// // We assume perm is allocated but not initialized
	// for(int i = 0; i < n_nodes; i++) 
	// 	perm[i] = i;

	// qsort_r(perm, n_nodes, sizeof(int), coord, cmp);

	// Create Node structs
	Node * nodes = malloc(n_nodes * sizeof(Node));
	for(int i = 0; i < n_nodes; i++) {
		nodes[i].i = i;
		nodes[i].x = coord[2*i];
		nodes[i].y = coord[2*i+1];
	}

	// Sort nodes
	qsort(nodes, n_nodes, sizeof(Node), cmp);

	// Fetch permutation (we assume perm is allocated)
	for(int i = 0; i < n_nodes; i++)
		perm[i] = nodes[i].i;

	return 0;
}