#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "romberg.h"

double romberg(double (*f)(double x), double a, double b, int n) {
    int kmax=1;
    double h = b-a; 
    double* R = malloc(n*n*sizeof(double));
    R[0] = (h/2.0)*(f(a)+f(b));

    for (int i=1; i<=n; i++) {

        h = h/2.0;
        double sum = 0;
        kmax = kmax*2;
        
        for (int k=1; k<=kmax-1; k+=2) {
            sum += f(a+k*h);
        }

        R[i*n] = 0.5 * R[(i-1)*n] + sum*h;
        
        for(int j=1; j<=i; j++){
            R[i*n+j] = R[i*n+j-1] + (R[i*n+j-1]-R[(i-1)*n+j-1]) / (pow(4.0, (double) j) - 1.0);
        }

    }
    
//    for (int l=0; l<n; l++){
//        for (int m=0; m<n; m++){
//            printf("%7.6f ", R[l*n+m]);
//        }
//        printf("\n");
//    }

    return R[n*n];
}

double f1(double x){
    return sin(sqrt(x));
}

double f2(double x){
    return sin(sqrt(x))-sqrt(x);
}
