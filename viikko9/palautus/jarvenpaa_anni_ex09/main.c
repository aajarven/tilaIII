#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"

int main(){
    double *xi = malloc(100*sizeof(double));
    double *xii1 = malloc(100*sizeof(double));
    double *xii2 = malloc(100*sizeof(double));
    double *yi = malloc(100*sizeof(double));
    double *yii1 = malloc(100*sizeof(double));
    double *yii2 = malloc(100*sizeof(double));

    srand(715517);

    for(int i=0; i<100; i++){
        xi[i] = i+1;
        yi[i] = xi[i]*0.5;
        xii1[i] = rand();
        xii2[i] = rand();
        yii1[i] = rand();
        yii2[i] = rand();
    }

    printf("2i)  %f\n", pearson(xi, yi, 100));
    printf("2ii-1) %f\n", pearson(xii1, yii1, 100));
    printf("2ii-2) %f\n", pearson(xii2, yii2, 100));
    printf("3i)  %f\n", kendall(xi, yi, 100));
    printf("3ii-1) %f\n", kendall(xii1, yii1, 100));
    printf("3ii-2) %f\n", kendall(xii2, yii2, 100));

//    double *sinX = malloc(2000*sizeof(double));
//    double *sinY = malloc(2000*sizeof(double));
//    for(int i=0; i<2000; i++){
//        sinX[i]=0.02*i;
//        sinY[i]=sin(sinX[i]);
//    }
//
//    printf("i\tautoc\n");
//    for(int i=0; i<100; i++){
//        printf("%d\t%f\n", i, autoc(xii, 100, i));
//    }

}
