#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "boundaryproblem.h"

#define PI (4*atan(1))

int main(){
    double h = PI/2.0/127.0;

    FILE *fp = fopen("shooting.dat", "w");
    shooting(0, PI/2.0, 1.0, 0, h, true, fp);


}
