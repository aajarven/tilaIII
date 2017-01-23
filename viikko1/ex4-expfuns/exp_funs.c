#include <stdio.h>
#include <math.h>

void exp_funs(double x){
    double a = (exp(x)-1)/x;
    double b = (exp(x)-exp(-x))/(2*x);

    printf("a) f1(%.6f)=%.6f\n", x, a)
    printf("b) f2(%.6f)=%.6f\n", x, b)
}
