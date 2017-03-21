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
    double* arr = fork_2d(0.5, 1.5, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.5, 1.5, arr[0], arr[1]);
    arr = fork_2d(0.4, 1.4, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", 0.4, 1.4, arr[0], arr[1]);
    arr = fork_2d(0.6, 1.6, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", -0.6, 1.6, arr[0], arr[1]);
    arr = fork_2d(-2, 4, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", -2.0, 4.0, arr[0], arr[1]);
    arr = fork_2d(-3, 5, &g);
    printf("[%f, %f]: (%.10f, %.10f)\n", -3.0, 5.0, arr[0], arr[1]);
}

double f(double x){ 
        return sqrt(x)+sin(x-sqrt(x));
}

double g(double x, double y){ 
        return 100*pow((y-x*x), 2)+pow(1-x, 2); 
}

