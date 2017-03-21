#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "integral.h"

double simpson(double a, double b, double eps, int level, int level_max){
    double h = b-a;
    double c = (a+b)*0.5;
    double one_simpson = h*(f(a)+4.0*f(c)+f(b))/6.0;
    double d = 0.5*(a+c);
    double e = 0.5*(c+b);
    double two_simpson = h*(f(a)+4.0*f(d)+2.0*f(c)+4.0*f(e)+f(b))/12.0;

    level++;
    if (level > level_max){
        return two_simpson;
    }

    if (fabs(two_simpson-one_simpson) < 15.0*eps) {
        return two_simpson + (two_simpson-one_simpson)/15.0;
    } else {
        double left_simpson = simpson(a, c, eps/2.0, level, level_max);
        double right_simpson = simpson(c, b, eps/2.0, level, level_max);
        return left_simpson + right_simpson;
    }
}

double f(double x){
    //return exp(-x)*cos(x)*cos(x);
    if (x==1){
        return 0;
    } else {
        return exp(x/(x-1))*cos(x/(1-x))*cos(x/(1-x))/((1-x)*(1-x));
    }
}
