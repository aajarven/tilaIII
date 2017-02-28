#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "re_der.h"

double re_der(double x){
   int N=5;
   double h=0.1;
   double D[N+1][N+1];

   for(int i=0; i<=N; i++){
       D[i][0]=(f(x+h)-f(x-h))/(2*h);
//       printf("(%d, %d)\n", i, 0);
       for(int j=0; j<=i-1; j++){
           D[i][j+1] = D[i][j] + (D[i][j]-D[i-1][j])/(pow(4.0,j+1.0)-1.0);
 //          printf("(%d, %d) (%d, %d) (%d, %d) \n", i, j+1, i, j, i-1, j);
       }
       h = h/2;
   }

//    for (int i=0; i<N; i++){
//        for(int j=0; j<=i; j++){
//            printf("%.10f\t", D[i][j]);
//        }
//        printf("\n");
//    }

   return D[N][N];
}

double f(double x){
    return sin(exp(x));
}
