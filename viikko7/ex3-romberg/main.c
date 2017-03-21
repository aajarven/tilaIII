#include <stdio.h>
#include <stdlib.h>
#include "romberg.h"

int main(){
    printf("A) n=3:  %.10f\n", romberg(&f1, 0, 1, 3));
    printf("   n=8:  %.10f\n", romberg(&f1, 0, 1, 8));
    printf("   n=10: %.10f\n", romberg(&f1, 0, 1, 10));
    printf("B) n=3:  %.10f\n", romberg(&f2, 0, 1, 3)+2.0/3.0);
    printf("   n=8:  %.10f\n", romberg(&f2, 0, 1, 8)+2.0/3.0);
    printf("   n=10: %.10f\n", romberg(&f2, 0, 1, 10)+2.0/3.0);
}
