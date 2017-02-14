#include <stdio.h>
#include <stdlib.h>
#include "roots.h"

int main(){
//    printf("bisect 0, 0.5: %.10f\n", bisect_f(0, 0.5));
//    printf("newton 0.25: %.10f\n", newton_f(0.25));
//    printf("bisect 0.5, 0,7: %.10f\n", bisect_f(0.5, 0.7));
//    printf("newton 0.6: %.10f\n", newton_f(0.6));

//    printf("newton 0, 0.1: %.10f\n", newton_g(0, 0.1));
//    printf("newton 0, 1: %.10f\n", newton_g(0, 1));
//    printf("newton 0, 10: %.10f\n", newton_g(0, 10));
//    printf("newton 0, 100: %.10f\n", newton_g(0, 100));

    double p[4] = {1, -6, 10, 3};
    double *roots = myroots(4, p);
    for (int i=0; i<3; i++){
        printf("%f + %f*i\n", roots[2*i], roots[2*i+1]);
    }

}
