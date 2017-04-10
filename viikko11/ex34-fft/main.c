#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fft.h"

#define PI (4*atan(1))

int main(){
    double h = 4.0*PI/127;

    //FILE *fp = fopen("fft.dat", "w");

    double *data = malloc((2*128+1)*sizeof(double));
    for (int i=1; i<2*128+1; i+=2){
        data[i] = sin(i*h);
        data[i+1] = 0;
    }

    fft(data, (unsigned int) 128, 1);

    for(int i=1; i<2*128+1; i+=2){
        printf("%f\t%f\n", data[i], data[i+1]);
    }


}
