#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double df(double, double);

double df(double x, double y){
    return x-y;
}

double euler(double x, double xf, double y, double h){
    while (xf-x > 0.00001){ // when testing x>fx step 0.1 ran one step too far
        double dy1 = df(x, y);
        double y2 = y + h*dy1;
        double dy2 = df(x+h, y2);
        
        y = y + 0.5*h*(dy1 + dy2);
        x += h;
    }

    return y;
}


double runge(double x, double xf, double y, double h){
    while(xf-x > 0.00001){
        double k1 = df(x, y);
        double k2 = df(x+h*0.5, y+h*0.5*k1);
        double k3 = df(x+h*0.5, y+h*0.5*k2);
        double k4 = df(x+h, y+h*k3);

        x += h;
        y += h/6.0*(k1+2*k2+2*k3+k4);
    }

    return y;
}

double adams3(double *x_in, double *y_in, double xf, double h){
    double x[3] = {x_in[0], x_in[1], x_in[2]};
    double y[3] = {y_in[0], y_in[1], y_in[2]};

    while(xf-x[0] > 0.0001){
        double x_next = x[0] + h;
        double y_next = y[0] + h/12.0*(23*df(x[0], y[0]) - 16*df(x[1], y[1]) +5*df(x[2], y[2]));
        y_next = y[0] + h/24.0*(9*df(x_next, y_next) + 19*df(x[0], y[0]) - 5*df(x[1], y[1]) + df(x[2], y[2]));

        y[2] = y[1];
        y[1] = y[0];
        y[0] = y_next;
        x[2] = x[1];
        x[1] = x[0];
        x[0] = x_next;

    }

    return y[0];
}

