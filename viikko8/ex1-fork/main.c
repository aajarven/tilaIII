#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forking.h"

double f(double);
double g(double, double);

int main(){
//    printf("[5.0, 10.0] %.10f\n", fork_1d(5.0, 10.0, &f));
//    printf("[5.1, 10.1] %.10f\n", fork_1d(5.1, 10.1, &f));
//    printf("[4.9, 9.9]  %.10f\n", fork_1d(4.9, 9.9, &f));
    double* mins = malloc(2*sizeof(double));
    double* maxes = malloc(2*sizeof(double));
    mins[0] = 0.0;
    mins[1] = 0.0;
    maxes[0] = 1.5;
    maxes[1] = 1.5;
    double* arr = fork_2d(mins, maxes, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.0, 1.5, arr[0], arr[1]);
    
    mins[0] = 0.4;
    mins[1] = 0.4;
    maxes[0] = 1.4;
    maxes[1] = 1.4;
    arr = fork_2d(mins, maxes, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.4, 1.4, arr[0], arr[1]);
    
    mins[0] = -0.6;
    mins[1] = -0.6;
    maxes[0] = 1.6;
    maxes[1] = 1.6;
    arr = fork_2d(mins, maxes, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", -0.6, 1.6, arr[0], arr[1]);

    mins[0] = 0.0;
    mins[1] = 0.0;
    maxes[2] = 10.0;
    maxes[1] = 10.0;
    arr = fork_2d(mins, maxes, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.0, 10.0, arr[0], arr[1]);
    
    mins[0] = 0.0;
    mins[1] = 0.0;
    maxes[0] = 20.0;
    maxes[1] = 20.0;
    arr = fork_2d(mins, maxes, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.0, 20.0, arr[0], arr[1]);
}

double f(double x){ 
        return sqrt(x)+sin(x-sqrt(x));
}

double g(double x, double y){ 
        return 100*pow((y-x*x), 2)+pow(1-x, 2); 
}

