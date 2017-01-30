#include <stdio.h>
#include <limits.h>
#include "harmonic_kahan.h"

int main(){
//    float summa = harmonic_kahan(1);
//    float edellinenSumma = 0;

    printf("%.7f\n", harmonic_kahan(INT_MAX-1));
//    printf("%.7f\n", harmonic_kahan(INT_MAX));

//    do {
//        edellinenSumma = summa;
//        summa = harmonic_kahan(N);
//        if(N%10000 == 0){
//            printf("%.6f\t%d\n", summa, N);
//        }
//        N+=1;
//    } while (edellinenSumma != summa);
//
//    printf("%.6f\n", summa);
}
