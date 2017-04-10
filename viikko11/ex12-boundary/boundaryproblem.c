#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define PI (4*atan(1))

double dy(double, double);
double runge(double, double, double, double, double, bool, FILE*);

double dy(double x, double y){
    return x-y;
}

/*
 * Calculates the shooting angle required to fulfill y''+y=x with boundary values
 * y(x0)=y0, y(xf)=yf using shooting method and RK4 with step h. Print defines whether
 * best found curve is printed.
 */
double shooting(double x0, double xf, double y0, double yf, double h, bool print, FILE *out){
    double k0 = -PI/6;
    double k1 = PI/6;

    double yk1 = runge(x0, xf, y0, k1, h, false, NULL);

    while(fabs(yk1-yf)>1e-10){
        double yk0 = runge(x0, xf, y0, k0, h, false, NULL);
        yk1 = runge(x0, xf, y0, k1, h, false, NULL);
        double k2 = k1 - (k1-k0)/(yk1-yk0)*yk1;
        k0 = k1;
        k1 = k2;
    }

    if (print){
        runge(x0, xf, y0, k1, h, true, out);
    }

    return k1;

}

double runge(double x, double xf, double y, double u, double h, bool print, FILE *out){
    while(xf-x > 1e-10){
        double k1y = h*u;
        double k1u = h*dy(x, y);
        double k2y = h*(u+.5*k1u);
        double k2u = h*dy(x+h*.5, y+.5*k1y);
        double k3y = h*(u+.5*k2u);
        double k3u = h*dy(x+h*.5, y+.5*k2y);
        double k4y = h*(u+k3u);
        double k4u = h*dy(x+h, y+k3y);

        x += h;
        y += 1/6.0*(k1y+2*k2y+2*k3y+k4y);
        u += 1/6.0*(k1u+2*k2u+2*k3u+k4u);

        if(print){
            fprintf(out, "%f\t%f\n", x, y);
        }
    }

    return y;
}

