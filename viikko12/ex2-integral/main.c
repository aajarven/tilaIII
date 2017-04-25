#include <float.h>
#include <stdio.h>
#include <time.h>
#include "integration.h"

int main(){
    double integral;
    clock_t start_t;
    clock_t end_t;

    printf("Rectangle rule:\nend\tsteps\tresult\ttime\n");

    for(double end=10; end<1e7; end=end*10){
        for(int steps=10; steps<1e9; steps=steps*10){
            start_t = clock();
            for(int i=0; i<10; i++){
                integral = rectRule(1.0, end, steps);
            }
            end_t = clock();
            printf("%e\t%e\t%f\t%e\n", end, (double)steps, integral, (double)(end_t-start_t)/CLOCKS_PER_SEC/10.0);
        }
    }


    printf("\n\nTrapezoidal rule:\nend\tsteps\tresult\ttime\n");

    for(double end=10; end<1e7; end=end*10){
        for(int steps=10; steps<1e9; steps=steps*10){
            start_t = clock();
            for(int i=0; i<10; i++){
                integral = trapRule(1.0, end, steps);
            }
            end_t = clock();
            printf("%e\t%e\t%f\t%e\n", end, (double)steps, integral, (double)(end_t-start_t)/CLOCKS_PER_SEC/10.0);
        }
    }


    printf("\n\nMonte Carlo:\nend\tsteps\tresult\ttime\n");

    for(double end=10; end<1e7; end=end*10){
        for(int steps=10; steps<1e9; steps=steps*10){
            start_t = clock();
            for(int i=0; i<10; i++){
                integral = monteCarlo(1.0, end, steps);
            }
            end_t = clock();
            printf("%e\t%e\t%f\t%e\n", end, (double)steps, integral, (double)(end_t-start_t)/CLOCKS_PER_SEC/10.0);
        }
    }
}
