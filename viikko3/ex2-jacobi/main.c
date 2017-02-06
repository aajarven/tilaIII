#include <stdio.h>
#include <stdlib.h>
#include "jacobi.h"
#include "lapack.h"

int main(){
    int size = 4;
    double mat[4*4] = {
        4, 30, 432, -543,
        30, 300, -675, 420,
        432, -675, 1620, -1050,
        -543, 420, -1050, 700
    };


    double *eigenvalues = jacobi(mat, size);
    printf("jacobi method:\n");
    for (int i=0; i<size; i++){
        printf("%f\t", eigenvalues[i]);
        if (i%size== size){
            printf("\n");
        }
    }

    char eigenvectors = 'N';
    double *WR = malloc(size*sizeof(double));
    double *WI = malloc(size*sizeof(double));
    
    // query best lwork
    int lwork = -1;
    double wkopt = 0;
    int info;
    dgeev_(&eigenvectors, &eigenvectors, &size, mat, &size, WR, WI, NULL, &size, NULL, &size, &wkopt, &lwork, &info);

    // actually run dgeev
    lwork = (int) wkopt;
    double *work = malloc(lwork*sizeof(double));
    dgeev_(&eigenvectors, &eigenvectors, &size, mat, &size, WR, WI, NULL, &size, NULL, &size, work, &lwork, &info);

    // print eigenvalues
    printf("\nLAPACK implementation:\n");
    for (int i=0; i<size; i++){
        printf("%f\t", WR[i]);
        if (i%size== size){
            printf("\n");
        }
    }
   
   printf("\n"); 
}
