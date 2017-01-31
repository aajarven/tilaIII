#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lapac.h"

double residual(int N, double *A, double *x, double *b, int m){

    double *C = (double*) malloc((size_t) N*sizeof(double));

    char TRANS = 'T';
    char NOTRANS = 'N';
    int K = N;
    double alpha = 1;
    double beta = 0;
    int M = N;
    int NN = 1;
    int LDA = N;
    int LDB = N;
    int LDC =N;

    int i;

    /* product Ax */
    dgemm_(&TRANS, &NOTRANS, &M, &NN, &K, &alpha, A, &LDA, x, &LDB, &beta, C, &LDC);
    
    /* Ax-b */
    for (i=0; i<N; i++){
        C[i] = C[i]-b[i];
    }

    /* residual */
    if (m==0){
        double biggest = fabs(C[0]);

        for (i=0; i<N; i++){
            biggest = fabs(C[i]) > biggest ? fabs(C[i]) : biggest;
        }

        return biggest;

    } else {
        double sum = 0;

        for (i=0; i<N; i++){
            sum += pow(fabs(C[i]), (double) m);
        }

        return pow(sum, 1.0/m);
    } 
}
