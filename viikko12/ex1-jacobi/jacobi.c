#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool doubleArrayEquals(int, double*, double*, double);

double* jacobi(int N, double *A, double *b, double *x0){
    double *x = malloc(N*sizeof(double));
    double *x_prev = malloc(N*sizeof(double));
    for (int i=0; i<N; i++){
        x[i] = x0[i];
    }

    do{
        for(int i=0; i<N; i++){
            x_prev[i] = x[i];
            printf("%e\t", x[i]);
        }
        printf("\n");

        for (int i=0; i<N; i++){
            double sigma = 0;
            for (int j=0; j<N; j++){
                if (j != i){
                    sigma += A[i*N+j]*x_prev[j];
                }
            }
            x[i] = 1.0/A[i*N+i]*(b[i]-sigma);
        }
    } while(!doubleArrayEquals(N, x_prev, x, 1e-12));

    return x;

}


bool doubleArrayEquals(int size, double *A1, double *A2, double tolerance){
    for(int i=0; i<size; i++){
        if (fabs(A1[i]-A2[i]) > tolerance){
            return false;
        }
    }

    return true;
}
