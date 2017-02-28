#include <stdio.h>
#include <stdlib.h>
#include "broyden.h"

int main(){
    double x_init[2] = {1.44, .14}; 
    double *x;
    for(int i=0; i<1000000; i++){
        x = broyden(x_init);
    }
    printf("x1 = %.10f\n", x[0]);
    printf("x2 = %.10f\n", x[1]);
}
