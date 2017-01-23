#include <stdio.h>
#include "harmonic.h"

int main(){
    printf("harmonic: %.6f\n", harmonic());

    printf("harmonic_bunch 50: %.6f\n", harmonic_bunch(50));
    printf("harmonic_bunch 100: %.6f\n", harmonic_bunch(100));
    printf("harmonic_bunch 500: %.6f\n", harmonic_bunch(500));
}
