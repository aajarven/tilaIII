#include <math.h>
#include <stdio.h>
#include "roots.h"

double bisect_f(double a, double b){
    double raja = 1e-10;
    int n=0;
    
    double alku = a;
    double loppu = b;
    double keskikohta = (a+b)/2.0;

    while (fabs(f(keskikohta)) > raja){
        if (signbit(f(alku)) == signbit(f(keskikohta))){
            alku = keskikohta;
        } else {
            loppu = keskikohta;
        }

        keskikohta = (alku+loppu)/2.0;
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return keskikohta;
}

double newton_f(double x0){
    int n=0;
    double raja = 1e-10;
    double x = x0;
    while (fabs(f(x))>raja){
        x = x-f(x)/df(x);
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return x;
}

double f(double x){
    double PI = acos(-1.0);
    return sin(3*PI*pow(x, 3)/(x*x-1))+.5;
}

double df(double x){
    double PI = acos(-1.0);
    return (9*PI*x*x/(x*x-1)-6*PI*pow(x, 4)/pow(x*x-1,2))*cos(3*PI*pow(x,3)/(x*x-1));
}

double newton_g(double x0, double B){
    int n=0;
    double raja = 1e-10;
    double x = x0;
    while (fabs(g(x, B))>raja){
        x = x-g(x, B)/dg(x, B);
        n++;
    }

    printf("result acquired after %d iterations\n", n);
    return x;
}

double g(double x, double B){
    return x+exp(-B*x*x)*cos(x);
}

double dg(double x, double B){
    return exp(-B*x*x)*(exp(B*x*x)-2*B*x*cos(x)-sin(x));
}
