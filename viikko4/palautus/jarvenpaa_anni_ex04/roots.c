#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "roots.h"
#include "lapack.h"

double bisect_f(double a, double b){
    double raja = 1e-10;
    int n=0;
    
    double alku = a;
    double loppu = b;
    double keskikohta = (a+b)/2.0;

    while (fabs(f(keskikohta)) > raja){
        if (signbit(f(alku)) == signbit(f(keskikohta))){
            alku = keskikohta;
        } else {
            loppu = keskikohta;
        }

        keskikohta = (alku+loppu)/2.0;
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return keskikohta;
}

double newton_f(double x0){
    int n=0;
    double raja = 1e-10;
    double x = x0;
    while (fabs(f(x))>raja){
        x = x-f(x)/df(x);
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return x;
}

double f(double x){
    double PI = acos(-1.0);
    return sin(3*PI*pow(x, 3)/(x*x-1))+.5;
}

double df(double x){
    double PI = acos(-1.0);
    return (9*PI*x*x/(x*x-1)-6*PI*pow(x, 4)/pow(x*x-1,2))*cos(3*PI*pow(x,3)/(x*x-1));
}

double newton_g(double x0, double B){
    int n=0;
    double raja = 1e-10;
    double x = x0;
    while (fabs(g(x, B))>raja){
        x = x-g(x, B)/dg(x, B);
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return x;
}

double g(double x, double B){
    return x+exp(-B*x*x)*cos(x);
}

double dg(double x, double B){
    return exp(-B*x*x)*(exp(B*x*x)-2*B*x*cos(x)-sin(x));
}

double* myroots(int N, double *p){
    int size = N-1;
    double *matriisi = malloc(size*size*sizeof(double));

    // first row
    for (int i=0; i<size; i++){
        matriisi[i] = -p[size-i-1]/p[size];
    }

    //other rows
    for (int i=0; i<size-1; i++){
       for (int j=0; j<size; j++){
           if (i==j){
               matriisi[(i+1)*size+j] = 1;
           } else {
               matriisi[(i+1)*size+j] = 0;
           }
        }
    }

//    for (int i=0; i<N-1; i++){
//        for (int j=0; j<N-1; j++){
//            printf("%f\t", matriisi[i*size+j]);
//        }
//        printf("\n");
//    }
    

    // eigenvector thingies
    char eigenvectors = 'N';
    double *WR = malloc((N-1)*sizeof(double));
    double *WI = malloc((N-1)*sizeof(double));
    int lwork = -1;
    double wkopt = 0;
    int info;
    
    // find best lwork
    dgeev_(&eigenvectors, &eigenvectors, &size, matriisi, &size, WR, WI, NULL, &size, NULL, &size, &wkopt, &lwork, &info);

    // actually find eigenvalues
    lwork = (int) wkopt;
    double *work = malloc(lwork*sizeof(double));
    dgeev_(&eigenvectors, &eigenvectors, &size, matriisi, &size, WR, WI, NULL, &size, NULL, &size, work, &lwork, &info);

    // construct return array
    double* ret = malloc(size*2*sizeof(double));
    for (int i=0; i<N-1; i++){
        ret[2*i] = WR[i];
        ret[2*i+1] = WI[i];
    }

    return ret;
}
