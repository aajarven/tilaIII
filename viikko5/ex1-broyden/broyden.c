#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "broyden.h"

double* broyden(double* x_init){
    // J
    double x1 = x_init[0];
    double x2 = x_init[1];
    double *J = malloc(4*sizeof(double));
    J[0] = -2*x1*exp(-x1*x1-x2*x2);
    J[1] = -2*x2*exp(-x1*x1-x2*x2);
    J[2] = cos(x1);
    J[3] = sin(x2);

    // B0 from inverse
    double *B = malloc(4*sizeof(double));
    double kerroin = 1/(J[0]*J[3]-J[1]*J[2]);

    B[0] = J[3]*kerroin;
    B[1] = -1*J[1]*kerroin;
    B[2] = -1*J[2]*kerroin;
    B[3] = J[0]*kerroin;
    
    // create arrays etc that are needed
    double *x = malloc(2*sizeof(double));
    x[0] = x1;
    x[1] = x2;
    double *x_next = malloc(2*sizeof(double));
    double *s = malloc(2*sizeof(double));
    double *y = malloc(2*sizeof(double));
    double *B_next = malloc(4*sizeof(double));
    double *f = malloc(2*sizeof(double));
    double *f_next = malloc(2*sizeof(double));
    double diff = 0;

    do{
        // 1
        f[0] = f1(x[0], x[1]);
        f[1] = f2(x[0], x[1]);
        double *tmpprod = matvecmul(B, f);

        //printf("B:\n%e\t%e\n%e\t%e\n\n", B[0], B[1], B[2], B[3]);
        //printf("f:\n%e\t%e\n\n", f[0], f[1]);
        //printf("tmpprod:\n%e\t%e\n\n", tmpprod[0], tmpprod[1]);
        x_next = vecdiff(x, tmpprod);
        //printf("x_next:\n%e\t%e\n\n", x_next[0], x_next[1]);
        
        // 2
        s = vecdiff(x_next, x);
        //printf("s:\n%e\t%e\n\n", s[0], s[1]);

        // 3
        f_next[0] = f1(x_next[0], x_next[1]);
        f_next[1] = f2(x_next[0], x_next[1]);
        //printf("f_next:\n%e\t%e\n\n", f_next[0], f_next[1]);
        
        y = vecdiff(f_next, f);
        //printf("y:\n%e\t%e\n\n", y[0], y[1]);

        // 4
        tmpprod = vecdiff(s, matvecmul(B, y));
        //printf("tmpprod:\n%e\t%e\n\n", tmpprod[0], tmpprod[1]);
        tmpprod = veccvecrmul(tmpprod, s);
        //printf("tmpprod:\n%e\t%e\n%e\t%e\n\n", tmpprod[0], tmpprod[1], tmpprod[2], tmpprod[3]);
        tmpprod = matmatmul(tmpprod, B);
        //printf("osoittaja:\n%e\t%e\n%e\t%e\n\n", tmpprod[0], tmpprod[1], tmpprod[2], tmpprod[3]);

        double denominator = vecrveccmul(vecmatmul(s, B), y);
        //printf("nimittäjä:\n%e\n\n", denominator);
        B_next[0] = B[0]+tmpprod[0]/denominator;
        B_next[1] = B[1]+tmpprod[1]/denominator;
        B_next[2] = B[2]+tmpprod[2]/denominator;
        B_next[3] = B[3]+tmpprod[3]/denominator;
        memcpy(B, B_next, 4*sizeof(double));
        diff = squaredifference(x, x_next);
        memcpy(x, x_next, 2*sizeof(double));

    } while (diff>1e-10);

    return x;

}

double f1(double x1, double x2){
    return exp(-(x1*x1+x2*x2))-0.125;
}

double f2(double x1, double x2){
    return sin(x1)-cos(x2);
}

/** 
 * v*M where M is 2*2 matrix and v is vector with length 2
 **/
double* vecmatmul(double *v, double *M){
    double *ret = malloc(2*sizeof(double));
    ret[0] = M[0]*v[0]+M[1]*v[0];
    ret[1] = M[2]*v[1]+M[3]*v[1];
    return ret;
}

/**
 * M*v where M is 2*2 and v is vector with length 2
 **/
double* matvecmul(double *M, double *v){
    double *ret = malloc(2*sizeof(double));
    ret[0] = M[0]*v[0] + M[1]*v[1];
    ret[1] = M[2]*v[0] + M[3]*v[1];
    return ret;
}

/**
 * v1-v2 where v1 and v2 are vectors with length 2
 * */
double* vecdiff(double *v1, double *v2){
    double *ret = malloc(2*sizeof(double));
    ret[0] = v1[0]-v2[0];
    ret[1] = v1[1]-v2[1];
    return ret;
}

/**
 * calculates sum of (M_1[i]-M_2[i])^2 when M_1 and M_2 are vectors with length 2
 * */
double squaredifference(double *M1, double *M2){
    double ret = 0;
    for (int i=0; i<2; i++){
        double diff = M1[i]-M2[i];
        ret += diff*diff;
    }
    return ret;
}

/**
 * multiplies column vector v1 by row vector v2
 * */
double* veccvecrmul(double *v1, double *v2){
    double *ret = malloc(4*sizeof(double));
    ret[0] = v1[0]*v2[0];
    ret[1] = v1[0]*v2[1];
    ret[2] = v1[1]*v2[0];
    ret[3] = v1[1]*v2[1];
    return ret;
}

/**
 * multiplies two 2*2 matrices
 * */
double* matmatmul(double *M1, double *M2){
    double *ret = malloc(4*sizeof(double));
    ret[0] = M1[0]*M2[0]+M1[1]*M2[2];
    ret[1] = M1[0]*M2[1]+M1[1]*M2[3];
    ret[2] = M1[2]*M2[0]+M1[3]*M2[2];
    ret[3] = M1[2]*M2[1]+M1[3]*M2[3];
    return ret;
}

/**
 * multiplies row vector v1 by column vector v2
 * */
double vecrveccmul(double *v1, double *v2){
    return v1[0]*v2[0]+v1[1]*v2[1];
}
