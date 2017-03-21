#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "forking.h"

// x and y to be used when moving only in direction of the other
double constX;
double constY;

// function given to fork_2d, stored in here because I need another function that only takes 1 parameter and calculates values of f(x, y) keeping either x or y constant. Ugliest hack possible but only thing I could make up without copypasting the code of fork_1d to fork_2d.
double (*f2)(double, double);


// note: naming conventions from lecture notes
double fork_1d(double a, double c, double (*f)(double)){
    double x;
    double e = sqrt(DBL_EPSILON);
    while (c-a > e){
        // using halfway points for b and x
        double b = (a+c)/2.0;
        x = (b+c)/2.0;
    
        if (f(x) < f(b)){
            a = b;
        } else {
            c = x;
        }
    }

    return x;
}

double* fork_2d(double a, double c, double (*f)(double, double)){
    // save given function to global
    f2 = f;

    // start from the middle
    constX = (a+c)/2;
    constY = (a+c)/2;

    // these are used to check wether we have converged
    double prevX = 0;
    double prevY = 0;
    
    double e = sqrt(DBL_EPSILON);
    while (abs(prevX-constX) > e || abs(prevY-constY) > e){
        // update previous values to current ones
        prevX = constX;
        prevY = constY;

        // calculate new ones
        constY = fork_1d(a, c, &fy);
        constX = fork_1d(a, c, &fx);
    }

    // return when converged
    double *ret = malloc(2*sizeof(double));
    ret[0] = constX;
    ret[1] = constY;
    return ret; 
}

double fy(double y){
    return f2(constX, y);
}

double fx(double x){
    return f2(x, constY);
}
