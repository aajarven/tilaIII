#include <stdio.h>
#include "harmonic_kahan.h"

int main(){
    float summa = harmonic_kahan(1);
    float edellinenSumma = 0;
    int N = 2;

    printf("%.6f\n", harmonic_kahan(10000000));
    printf("%.6f\n", harmonic_kahan(10000001));

    do {
        edellinenSumma = summa;
        summa = harmonic_kahan(N);
//        if(N%100 == 0){
//            printf("%.6f\n", summa);
//        }
        N+=1;
    } while (edellinenSumma != summa);

    printf("%.6f\n", summa);
}
