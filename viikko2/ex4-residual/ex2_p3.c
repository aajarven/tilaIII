#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  int n,i,j,c1,c2,*pivot,ok;
  double *A,*b;
  
  scanf("%d",&n);
  printf("\nn %d\n",n);
  A=(double *)malloc((size_t)n*n*sizeof(double));
  b=(double *)malloc((size_t)n*sizeof(double));
  pivot=(int *)malloc((size_t)n*sizeof(int));

  for (i=0;i<n;i++) for (j=0;j<n;j++) scanf("%lg",&A[j*n+i]);
  for (i=0;i<n;i++) scanf("%lg",&b[+i]);

  printf("A\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) printf("%12.8g ",A[j*n+i]);
    printf("\n");
    }
  printf("\nb\n");
  for (i=0;i<n;i++) printf("%12.8g ",b[i]);
  printf("\n");

  c1=n;
  c2=1;
  
  dgesv_(&c1, &c2, A, &c1, pivot, b, &c1, &ok);      
  
  printf("x\n");
  for (i=0;i<n;i++) printf("%12.8g ", b[i]);	
  printf("\n\n");

  return 0;
} 
