#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "re_der.h"

double erotusosamaara(int);

int main(){
//    double realDerivative = exp(1.0)*cos(exp(1.0));
//
    printf("re_der %.20f\n\n", re_der_h(1, 0.1, 5));
//    printf("re_der %.20f\n\n", re_der_h(1.0, 0.5, 6));
//
//    printf("n=2 %.20f - %.20e = %.20e\n", realDerivative, re_der_h(1, pow(10.0, -2)), (realDerivative-re_der_h(1, pow(10.0, -2.0))));
//
//
//    printf("\nabsoluuttinen virhe\n");
//    printf("n=2 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -2.0), 5)));
//    printf("n=4 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -4.0), 5)));
//    printf("n=6 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -6.0), 5)));
//    printf("n=8 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -8.0), 5)));
//    printf("n=10 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -10.0), 5)));
//    printf("n=12 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -12.0), 5)));
//    printf("n=14 %.10e\n", (realDerivative-re_der_h(1, pow(10.0, -14.0), 5)));
//
//    printf("\nsuhteellinen virhe\n");
//    printf("n=2 %.10e\n", (re_der_h(1, pow(10.0, -2.0), 5) - realDerivative) / realDerivative);
//    printf("n=4 %.10e\n", (re_der_h(1, pow(10.0, -4.0), 5) - realDerivative) / realDerivative);
//    printf("n=6 %.10e\n", (re_der_h(1, pow(10.0, -6.0), 5) - realDerivative) / realDerivative);
//    printf("n=8 %.10e\n", (re_der_h(1, pow(10.0, -8.0), 5) - realDerivative) / realDerivative);
//    printf("n=10 %.10e\n", (re_der_h(1, pow(10.0, -10.0), 5) - realDerivative) / realDerivative);
//    printf("n=12 %.10e\n", (re_der_h(1, pow(10.0, -12.0), 5) - realDerivative) / realDerivative);
//    printf("n=14 %.10e\n", (re_der_h(1, pow(10.0, -14.0), 5) - realDerivative) / realDerivative);
//
//    printf("\nerotusosamäärät\n");
//    for (int i=4; i<=14; i+=2){
//        printf("n=%d %.10e\n", i, erotusosamaara(i));
//    }
//    
//    printf("\nsuhteellinen virhe\n");
//    for (int i=4; i<=14; i+=2){
//        printf("n=%d %.10e\n", i, (erotusosamaara(i)-realDerivative)/realDerivative);
//    }
}

double erotusosamaara(int n){
    return (sin(exp(1+pow(10.0, -n)))-sin(exp(1-pow(10.0, -n))))/(2*pow(10.0, -n));
}
