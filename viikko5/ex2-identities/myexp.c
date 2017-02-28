#include<stdio.h>
#include<math.h>

double myexp(double x){
    double ln2= log(2);
    return pow(2,(round(x/ln2)))*pow(2,(x/ln2-round(x/ln2)));
}
