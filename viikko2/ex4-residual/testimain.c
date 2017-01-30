#include <stdio.h>
#include <stdlib.h>
//#include "harmonic_kahan.h"

int main(){
	double *A = (double*) malloc((size_t) 3*4*sizeof(double));
	double *B = (double*) malloc((size_t) 4*2*sizeof(double));
	double *C = (double*) malloc((size_t) 3*2*sizeof(double));
	
	int i;
	for (i=0; i<3*4; i++){
		A[i] = i;
	}

	for (i=0; i<4*2; i++){
		B[i] = 8-i;
	}

	char TRANS = 'T';
	int K = 4;
	double alpha = 1;
	double beta = 1;

	dgemm_(&TRANS, &TRANS, &K, &alpha, &beta, &C);
	
	int j;
	for (i=0; i<3; i++){
		for (j=0; j<2; j++){
			printf("%d\t", C[i*3+j]);
		}
	}

}
