#include <stdio.h>
#include <stdlib.h>
#include "integral.h"

int main(){
    printf("%.10f\n", simpson(0, 1, 1e-10, 0, 80));
}
