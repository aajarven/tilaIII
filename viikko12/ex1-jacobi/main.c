#include <stdio.h>
#include <time.h>
#include "jacobi.h"

int main(){
    double A[16] = {1.0, 2.0, 3.0, 4.0, 1.0, 3.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 3.0, 1.0, 0.0, 1.0};
    double b[4] = {20.0, 11.0, 6.0, 4.0};
    double x0[4] = {0.0+1e-19, 1.0+1e-19, 2.0+1e-19, 3.0+1e-19};

    double w = 16;

    for(int i=0; i<4; i++){
        A[i*4+i] += A[i*4+i]+w;
    }
    
    double *x = jacobi(4, A, b, x0);
    

    
    printf("x: [");
    for(int i=0; i<4; i++){
        printf("%f\t", x[i]);
    }
    printf("]\n");
}
