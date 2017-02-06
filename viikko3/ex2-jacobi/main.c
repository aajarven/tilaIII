#include <stdio.h>
#include <stdlib.h>
#include "jacobi.h"
#include "lapack.h"

int main(){
    int size = 4;
    double *mat = malloc(size*size*sizeof(double));
    mat[0] = 4;
    mat[1] = -30;
    mat[2] = 60;
    mat[3] = -35;
    mat[4] = -30;
    mat[5] = 300;
    mat[6] = -675;
    mat[7] = 420;
    mat[8] = 60;
    mat[9] = -675;
    mat[10]= 1620;
    mat[11]= -1050;
    mat[12]= -35;
    mat[13]= 420;
    mat[14]= -1050;
    mat[15]= 700;


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
