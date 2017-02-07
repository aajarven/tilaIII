#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "jacobi.h"

struct matrixIndex{
    int i;
    int j;
};


/**
 * Jacobi method as presented in http://www.cmi.ac.in/~ksutar/NLA2013/iterativemethods.pdf
 **/

double* jacobi(double *Q, int N){
    double PI = acos(-1.0);
    double threshold = 1e-10; // how close to zero should non-diagonal elements be

    // copy of Q
    double *D = malloc(N*N*sizeof(double));
    memcpy(D, Q, N*N*sizeof(double));

    // unit matrix
    double *S = malloc(N*N*sizeof(double));
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            if (i==j){
                S[i*N+j] = 1;
            } else {
                S[i*N+j] = 0;
            }
        }
    }

    // S_1 & S_1t
    double *S1 = malloc(N*N*sizeof(double));
    double *S1t = malloc(N*N*sizeof(double));

    bool diagonal = false;
    do{
        // find biggest off diagonal value
        struct matrixIndex biggest = biggestOffDiag(D, N);

        // finding the rotational angle
        double theta;
        if (D[biggest.i*N+biggest.i] == D[biggest.j*N+biggest.j]){
            if (D[biggest.i*N+biggest.j] > 0){
                theta = PI/4;
            } else {
                theta = -PI/4;
            }
        } else {
            theta = 0.5*atan(2*D[biggest.i*N+biggest.j]/(D[biggest.i*N+biggest.i]-D[biggest.j*N+biggest.j]));
        }

        // compute S_1 (and its transpose)
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                if (i==j){
                    S1[i*N+j] = 1;
                    S1t[i*N+j] = 1;
                } else {
                    S1[i*N+j] = 0;
                    S1t[i*N+j] = 0;
                }
            }
        }
        S1[biggest.i*N+biggest.i] = cos(theta);
        S1[biggest.j*N+biggest.j] = S1[biggest.i*N+biggest.i];
        S1[biggest.j*N+biggest.i] = sin(theta);
        S1[biggest.i*N+biggest.j] = -S1[biggest.j*N+biggest.i];
        S1t[biggest.i*N+biggest.i]=S1[biggest.i*N+biggest.i];
        S1t[biggest.j*N+biggest.j]=S1[biggest.j*N+biggest.j];
        S1t[biggest.i*N+biggest.j]=S1[biggest.j*N+biggest.i];
        S1t[biggest.j*N+biggest.i]=S1[biggest.i*N+biggest.j];

        // matrix for doing multiplication S1t*D*S1
        double *tmp = malloc(N*N*sizeof(double));

        // S1t*D
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                tmp[i*N+j] = 0;
                for (int k=0; k<N; k++){
                    tmp[i*N+j] += S1t[i*N+k]*D[k*N+j];
                }
            }
        }

        // *S1
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                D[i*N+j] = 0;
                for (int k=0; k<N; k++){
                    D[i*N+j] += tmp[i*N+k]*S1[k*N+j];
                }
            }
        }

        // S*S1
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                tmp[i*N+j] = 0;
                for (int k=0; k<N; k++){
                    tmp[i*N+j] += S[i*N+j]*S1[i*N+j];
                }
            }
        }

        // check diagonality
        diagonal = true;
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                if (i != j && fabs(D[i*N+j]) > threshold){
                    diagonal = false;
                }
            }
        }
    
        free(tmp);

    } while (!diagonal);

    double *eigenvalues = malloc(N*sizeof(double));
    for (int i=0; i<N; i++){
        eigenvalues[i] = D[i*N+i];
    }

    free(S);
    free(S1);
    free(S1t);
    free(D);

    return eigenvalues;

}

double err_propag(int N, double dq){
    double *M = malloc(N*N*sizeof(double));

    // fill with random numbers
    srand(time(0));
    for (int i=0; i<N; i++){
        for (int j=i; j<N; j++){
            M[i*N+j] = rand()/RAND_MAX;
            M[j*N+i] = M[i*N+j];
        }
    }

    double *eigen1 = jacobi(M, N);

    //perturb and recalculate
    int x = rand() % (N*N);
    M[x] = M[x]*dq;

    double *eigen2 = jacobi(M, N);

    double eigennorm = 0;
    double errornorm = 0;
    for (int i=0; i<N; i++){
        eigennorm += eigen1[i] * eigen1[i];
        errornorm += (eigen2[i]-eigen1[i]);
    }

    return pow(errornorm, 0.5)/(pow(eigennorm, 0.5)*dq);
}


struct matrixIndex biggestOffDiag(double *matrix, int size){

    struct matrixIndex ret = {0, 1};

    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            if (i != j && fabs(matrix[i*size+j]) > fabs(matrix[ret.i*size+ret.j])){
                ret.i = i;
                ret.j = j;
            }
        }
    }

    return ret;
}
