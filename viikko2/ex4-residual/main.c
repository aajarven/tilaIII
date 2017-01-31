#include <stdio.h>
#include <stdlib.h>
#include "residual.h"
#include "lapac.h"

int main(){

	int N, i, j, M;
	double *A, *B, *C;

        /* read N and M */
	scanf("%d",&N);
	scanf("%d",&M);

	A = (double*) malloc((size_t) N*N*sizeof(double));
	B = (double*) malloc((size_t) 1*N*sizeof(double));
	C = (double*) malloc((size_t) 1*N*sizeof(double));

	/* reading A */
	for (i=0;i<N;i++) for (j=0;j<N;j++) scanf("%lg",&A[j*N+i]);

	/* reading C (corresponds to b */
	for (i=0;i<N;i++) scanf("%lg",&C[+i]);
	
        /* reading B (corresponds to x */
	for (i=0;i<N;i++) scanf("%lg",&B[+i]);


	printf("%f\n", residual(N, A, B, C, M));
	return 1;
}
