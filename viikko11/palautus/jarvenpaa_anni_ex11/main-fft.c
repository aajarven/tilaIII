#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fft.h"

#define PI (4*atan(1))

int main(){
    double h = 4.0*PI/127;
    unsigned int n = 128;
    //FILE *fp = fopen("fft.dat", "w");

    double *data = malloc((2*n+1)*sizeof(double));
    srand(715517);
    for (unsigned int i=1; i<2*n+1; i+=2){
        // sin
        //data[i] = sin(i*h);
        
        // 3b
        //data[i] = sin(i*h) + 2.0*rand()/(double)(RAND_MAX)-0.5;
        
        data[i] = sin(10*i*h)+sin(10.5*i*h);
        data[i+1] = 0;
    }

    fft(data, n, 1);

    // truncate and scale data
    double *truncatedData = malloc(n*sizeof(double));
    for (unsigned int i=0; i<n/2; i++){
        truncatedData[i] = 2.0*data[2*i+1]/n;
        truncatedData[i+1] = 2.0*data[2*i+2]/n;
    }

    double *magnitudes = malloc(n/2*sizeof(double));
    for(unsigned int i=0; i<n/2; i++){
        magnitudes[i] = sqrt(pow(truncatedData[2*i], 2.0) + pow(truncatedData[2*i+1], 2.0));
    }

    double *frequencies = malloc(n/2*sizeof(double));
    for(unsigned int i=0; i<n/2; i++){
        frequencies[i] = i*h;
    }

    for(unsigned int i=0; i<n/2; i++){
        printf("%f\t%f\n", frequencies[i], magnitudes[i]);
    }




//    for(int i=1; i<2*128+1; i+=2){
//        printf("%f\t%f\n", data[i], data[i+1]);
//    }


}
