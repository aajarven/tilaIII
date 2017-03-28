#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"

int main(){
    double *xi = malloc(100*sizeof(double));
    double *xii = malloc(100*sizeof(double));
    double *yi = malloc(100*sizeof(double));
    double *yii = malloc(100*sizeof(double));

    srand(715517);

    for(int i=0; i<100; i++){
        xi[i] = i+1;
        yi[i] = xi[i]*0.5;
        xii[i] = rand();
        yii[i] = rand();
    }

    printf("2i)  %f\n", pearson(xi, yi, 100));
    printf("2ii) %f\n", pearson(xii, yii, 100));
    printf("3i)  %f\n", kendall(xi, yi, 100));
    printf("3ii) %f\n", kendall(xii, yii, 100));

    double *sinX = malloc(2000*sizeof(double));
    double *sinY = malloc(2000*sizeof(double));
    for(int i=0; i<2000; i++){
        sinX[i]=0.02*i;
        sinY[i]=sin(sinX[i]);
    }

//    printf("i\tautoc\n");
//    for(int i=0; i<100; i++){
//        printf("%d\t%f\n", i, autoc(xii, 100, i));
//    }

}
