#include <stdlib.h>


double f(double x){
    return 1.0/(x*x);
}

double rectRule(double start, double end, int steps){
    double h = (end-start)/steps;
    
    double integral = 0;
    for (double mid = start+h/2; mid<end; mid+=h){
        integral += f(mid)*h;
    }

    return integral;
}


double trapRule(double start, double end, int steps){
    double h = (end-start)/steps;

    double integral = 0;
    for(double substart = start; substart<end; substart+=h){
        integral +=(f(substart)+f(substart+h))/2.0*h;
    }

    return integral;
}


double monteCarlo(double start, double end, int steps){
    srand(715517);
    
    double integral = 0;
    double r;
    double length = end-start;
    for (int i=0; i<steps; i++){
        r = (double)rand()/(double)(RAND_MAX)*length+start;
        integral += f(r);
    }

    return integral*length/steps;
}
