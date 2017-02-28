#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "re_der.h"

double re_der(double x){
   int N=5;
   double h=0.1;
   double *D = malloc((N+1)*(N+1)*sizeof(double));

   for(int i=0; i<=N; i++){
       D[i*(N+1)]=(f(x+h)-f(x-h))/(2*h);
//       printf("%.6f\t", D[i*N]);
       for(int j=0; j<i-1; j++){
//           printf("(%d, %d), (%d, %d), (%d, %d)\n", i, j+1, i, j, (i-1), j);
//           printf("(%d, %d)\n", i, j);
           D[i*(N+1)+j+1] = D[i*(N+1)+j] + (D[i*(N+1)+j]-D[(i-1)*(N+1)+j])/(pow(4.0,j+1.0)-1.0);
//           printf("%.6f\t", D[i*N+j+1]);
       }
//       printf("\n");
       h = h/2;
   }

    for (int i=0; i<N; i++){
        for(int j=0; j<i; j++){
            printf("%.10f\t", D[i*(N+1)+j]);
        }
        printf("\n");
    }

   return D[(N+1)*(N+1)];
}

double f(double x){
    return sin(exp(x));
}
