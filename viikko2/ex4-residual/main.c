#include <stdio.h>
#include <stdlib.h>
#include "residual.h"
#include "lapac.h"

int main(){
	double *A = (double*) malloc((size_t) 4*4*sizeof(double));
	double *B = (double*) malloc((size_t) 1*4*sizeof(double));
	double *C = (double*) malloc((size_t) 1*4*sizeof(double));

	int i;
	for (i=0; i<4*4; i++){
		A[i] = i;
	}

	for (i=0; i<4; i++){
		B[i] = 8.0-i;
	}

	
	for (i=0; i<4; i++){
		C[i] = 1.0+i;
	}

	printf("\n%f\n", residual(4, A, B, C, 3));
	
	return 1;
}
