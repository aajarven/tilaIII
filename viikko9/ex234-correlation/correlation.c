#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "correlation.h"

double pearson(double* x, double* y, int N){
    double meanX = mean(x, N);
    double meanY = mean(y, N);

    double numerator = 0;
    double denominator1 = 0;
    double denominator2 = 0;
    for(int i=0; i<N; i++){
        numerator += (x[i]-meanX)*(y[i]-meanY);
        denominator1 += pow(x[i]-meanX, 2);
        denominator2 += pow(y[i]-meanY, 2);
    }

    return numerator/(sqrt(denominator1)*sqrt(denominator2));
}

double mean(double *arr, int len){
    double sum = 0;

    for(int i=0; i<len; i++){
        sum += arr[i];
    }

    return sum/len;
}

double kendall(double* x, double* y, int N){
    int concordant = 0;
    int discordant = 0;

    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if(x[i]<x[j] && y[i]<y[j]){
                concordant++;
            } else if(x[i]>x[j] && y[i]>y[j]){
                concordant++;
            } else {
                discordant++;
            }
        }
    }

    return (concordant - discordant)/(N*(N-1)/2.0);
}

double autoc(double* x, int N, int k){
    double meanX = mean(x, N);

    double numerator = 0;
    for(int i=k; i<N; i++){
        numerator += (x[i]-meanX)*(x[i-k]-meanX);
    }
    
    double denominator = 0;
    for(int i=0; i<N; i++){
        denominator += pow(x[i]-meanX, 2.0);
    }

    return numerator/denominator;
}
