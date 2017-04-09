#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double df(double, double);

double df(double x, double y){
    return x-y;
}

double shooting(double x0, double xf, double y, double h){
    while (xf-x > 0.00001){ // when testing x>fx step 0.1 ran one step too far
        double dy1 = df(x, y);
        double y2 = y + h*dy1;
        double dy2 = df(x+h, y2);
        
        y = y + 0.5*h*(dy1 + dy2);
        x += h;
    }

    return y;
}
