#include <stdio.h>
#include <time.h>
#include "odes.h"


int main(){
    double h[4] = {0.2, 0.1, 0.01, 0.001};
    
    for(int i=0; i<4; i++){
//        double result;
//        clock_t begin = clock();
//        for(int j=0; j<1000000; j++){
//            result = euler(0, 1, 0, h[i]);
//        }
//        clock_t end = clock();
        double result = euler(0, 1, 0, h[i]);
        printf("2a) h=%f:\t%f\n", h[i], result);
//        printf("    time spent on 1 000 000 runs:\t%f s\n", (double)(end-begin)/CLOCKS_PER_SEC);
    }
    
    for(int i=0; i<4; i++){
//        double result;
//        clock_t begin = clock();
//        for(int j=0; j<1000000; j++){
//            result = runge(0, 1, 0, h[i]);
//        }
//        clock_t end = clock();
        double result = runge(0, 1, 0, h[i]);
        printf("3) h=%f:\t%f\n", h[i], result);
//        printf("    time spent on 1 000 000 runs:\t%f s\n", (double)(end-begin)/CLOCKS_PER_SEC);
    }
}
