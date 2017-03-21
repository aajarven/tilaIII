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
        //printf("x: %f\n", x);
    }

    return x;
}

double* fork_2d(double *a, double *c, double (*f)(double, double)){
    // save given function to global
    f2 = f;

    // start from the middle
    constX = (a[0]+c[0])/2;
    constY = (a[1]+c[1])/2;

    //printf("constX: %f\nconstY: %f\n", constX, constY);
    
    // these are used to check wether we have converged
    double prevX = DBL_MAX;
    double prevY = DBL_MAX;
    //printf("prevX: %f\nprevY: %f\n", prevX, prevY);
    //printf("diffX: %f\ndiffY: %f\n", abs(prevX-constX), abs(prevY-constY));
    
    double e = sqrt(DBL_EPSILON);
    while (fabs(prevX-constX) > e || fabs(prevY-constY) > e){
        // update previous values to current ones
        prevX = constX;
        prevY = constY;

        //printf("before constX: %f\nconstY: %f\n", constX, constY);
        // calculate new ones
        constX = fork_1d(a[0], c[0], &fx);
        constY = fork_1d(a[1], c[1], &fy);
        //printf("constX: %f\nconstY: %f\n", constX, constY);
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
