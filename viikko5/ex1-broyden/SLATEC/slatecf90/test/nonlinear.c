
/* 
   Demonstration of calling the SLATEC routine DNSQE from C.
   
   Find the root of nonlinear equations:
      f1(x1,x2)=0, f2(x1,x2)=0
   where
      f1(x1,x2)=exp(-x1**2-x2**2)-1/8
      f2(x1,x2)=sin(x1)-cos(x2)

   Compilation:
      gcc nonlinear.c -lslatec -L../src -lgfortran
      gfortran nonlinear.c -lslatec -L../src

          Antti Kuronen, 2011
          2017: Corrected a bug: Jacobian was transposed
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 2
#define LWA ((3*N*N+13*N)/2)

void fcn(int *n, double *x, double *fvec, int *iflag) {
  fvec[0]=exp(-x[0]*x[0]-x[1]*x[1])-0.125;
  fvec[1]=sin(x[0])-cos(x[1]);
  if (*iflag==0) printf("%20.10g%20.10g%20.10g%20.10g\n",x[0],x[1],fvec[0],fvec[1]);
}

void jac(int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag) {
  /* 
     Now, be careful here and note the order of array fjac.  Fortran
     version has a 2D array and the order of assignments is
     11,12,21,22. However, Fortran memory layout convention says that
     the array is stored in the memory in order 11,21,12,22.
  */
  fjac[0]=-2.0*x[0]*exp(-x[0]*x[0]-x[1]*x[1]);
  fjac[1]=cos(x[0]);
  fjac[2]=-2.0*x[1]*exp(-x[0]*x[0]-x[1]*x[1]);
  fjac[3]=sin(x[1]);
}

int main (int argc, char **argv)
{
  double tol,x[N],fvec[N],wa[LWA];
  int n,lwa,iopt,nprint,info;

  n=N;
  lwa=LWA;

  x[0]=atof(*++argv);
  x[1]=atof(*++argv);
  tol=atof(*++argv);
  nprint=atoi(*++argv);
  iopt=1;

  dnsqe_(fcn,jac,&iopt,&n,x,fvec,&tol,&nprint,&info,wa,&lwa);
  printf("# finished, info=%d\n",info);
  printf("%20.10g%20.10g%20.10g%20.10g\n",x[0],x[1],fvec[0],fvec[1]);

  return(0);
}
