#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lu.h"

 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

 #define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


int lu(Matrix * A) {
	int m = A->m;
	Matrix * LU = A; // factorization is done in-place
	for(int k = 0; k < m-1; k++) {
		if(fabs(LU->a[k][k]) < EPS) return -1;
		for(int j = k+1; j < m; j++) {
			LU->a[j][k] /= LU->a[k][k];
			for(int l = k+1; l < m; l++)
				LU->a[j][l] -= LU->a[j][k] * LU->a[k][l];
		}
	}
	return 0;
}

int lu_band(BandMatrix * A) {
	int m = A->m, b = A->k;
	BandMatrix * LU = A; // factorization is done in-place
	for(int k = 0; k < m-1; k++) {
		if(fabs(LU->a[k][k]) < EPS) return -1;
		for(int j = k+1; j <= min(m-1, k+b); j++) {
			LU->a[j][k] /= LU->a[k][k];
			for(int l = k+1; l <= min(m-1, k+b); l++)
				LU->a[j][l] -= LU->a[j][k] * LU->a[k][l];
		}
	}
	return 0;
}

int solve(Matrix * LU, double * y) {
	int m = LU->m;
	double * x = y;

	// Résolution de L*x = y par substitution avant
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < i; j++)
			x[i] -= LU->a[i][j] * x[j];
	}

	// Résolution de U*x = L^{-1} y par substitution arrière
	for(int i = m-1; i >= 0; i--) {
		for(int j = i+1; j < m; j++)
			x[i] -= LU->a[i][j] * x[j];
		x[i] /= LU->a[i][i];
	}
	
	return 0;
}

int solve_band(BandMatrix * LU, double * y) {
	int m = LU->m, b = LU->k;
	double * x = y;

	// Résolution de L*x = y par substitution avant
	for(int i = 0; i < m; i++) {
		for(int j = max(0, i-b); j < i; j++) {
			x[i] -= LU->a[i][j] * x[j];
		}
	}

	// Résolution de U*x = L^{-1} y par substitution arrière
	for(int i = m-1; i >= 0; i--) {
		for(int j = i+1; j <= min(m-1, i+b); j++)
			x[i] -= LU->a[i][j] * x[j];
		x[i] /= LU->a[i][i];
	}

	return 0;
}