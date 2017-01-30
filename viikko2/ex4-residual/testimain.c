#include <stdio.h>
#include <stdlib.h>
/*#include "harmonic_kahan.h"*/

int main(){
	double *A = (double*) malloc((size_t) 3*4*sizeof(double));
	double *B = (double*) malloc((size_t) 4*2*sizeof(double));
	double *C = (double*) malloc((size_t) 3*2*sizeof(double));
	
	int i;
	for (i=0; i<3*4; i++){
		A[i] = i;
	}

    printf("A:\n");
	int j;
	for (i=0; i<3; i++){
		for (j=0; j<4; j++){
			printf("%f\t", A[i*4+j]);
		}
        printf("\n");
	}
	
    for (i=0; i<4*2; i++){
		B[i] = 8-i;
	}

    printf("\nB:\n");
	for (i=0; i<4; i++){
		for (j=0; j<2; j++){
			printf("%f\t", B[i*2+j]);
		}
        printf("\n");
	}
	
    char TRANS = 'C';
	int K = 4;
	double alpha = 1;
	double beta = 1;

	int M = 3;
    int N = 2;
    int LDA = 4;
    int LDB = 2;
    int LDC = 3;

	dgemm_(&TRANS, &TRANS, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
    
    printf("\nC:\n");

    // HUOM OBS NB sylkee transpoosi
	/*int j;*/
	for (i=0; i<2; i++){
		for (j=0; j<3; j++){
			printf("%f\t", C[i*3+j]);
		}
        printf("\n");
	}

}
