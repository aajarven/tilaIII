#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double df(double, double);

double df(double x, double y){
    return x-y;
}

double euler(double x, double xf, double y, double h){
    while (x < xf){
        double dy = df(x, y);
        x = x+h;
        y = y+h*dy;
    }

    return y;
}


double runge(double x, double xf, double y, double h){
    while(x < xf){
        double k1 = df(x, y);
        double k2 = df(x+h*0.5, y+h*0.5*k1);
        double k3 = df(x+h*0.5, y+h*0.5*k2);
        double k4 = df(x+h, y+h*k3);

        x += h;
        y += h/6.0*(k1+2*k2+2*k3+k4);
    }

    return y;
}
