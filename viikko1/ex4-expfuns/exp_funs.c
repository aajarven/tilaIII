#include <stdio.h>
#include <math.h>

void exp_funs(double x){
    double a = (exp(x)-1)/x;
    double b = (exp(x)-exp(-x))/(2*x);

     printf("%.6f\t%.6f\t%.6f\n", x, a, b);
//    printf("a) f1(%.6f)=%.6f\n", x, a);
//    printf("b) f2(%.6f)=%.6f\n", x, b);
}

void exp_funs_Taylor(double x){
    double a = 0;
    double prev_a;
    double i = 1;
    double factorial = 1; // current denominator
    do{
        prev_a = a;
        factorial *= i;
        a += pow(x, i-1)/factorial;
        i++;
    } while (a!=prev_a);

    a = a;

    double b = 0;
    double prev_b = 1;
    i = 1;
    factorial = 1;
    do{
        prev_b = b;
        factorial *= (i+1)*(i+2);
        b += pow(x, i-1)/factorial;
        i += 2;
    } while (b != prev_b);

    b = b;

     printf("%.6f\t%.6f\t%.6f\n", x, a, b);
//    printf("a) f1(%.6f)=%.6f\n", x, a);
//    printf("b) f2(%.6f)=%.6f\n", x, b);
}
