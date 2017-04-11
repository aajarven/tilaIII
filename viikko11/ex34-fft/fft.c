#include <math.h>
#include <stdio.h>

#define PI (4*atan(1))
#define SWAP(a,b) tmp=(a); (a)=(b); (b)=tmp

/*
 * Edits given array to contain discrete Fourier transform of given array if isign
 * equals 1 or if isign is -1 its inverse discrete Fourier transform.
 * Given data can be a complex array of length nn or a real array of length 2*nn.
 * In both cases nn must be an integer power of 2.
 *
 * Implemented following Numerical Recipes in C using Danielson-Lanczos
 */
void fft(double *data, unsigned int nn, int isign){
    double tmp, wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    unsigned int istep, n, mmax, m, i, j;

    n = nn<<1;
    j=1;
    for(i=1; i<n; i+=2){
        if (j>i){
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }

        m=nn;
        while(m>=2 && j>m){
            j -= m;
            m >>= 1;
        }

        j+=m;
    }

    mmax = 2;
    while (n > mmax){
        istep = mmax << 1;
        theta = isign*(2.0*PI/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;

        for(m=1; m<mmax; m+=2){
            for(i=m; i<=n; i+=istep){
                j = i+mmax;
                tempr = wr*data[j]-wi*data[j+1];
                tempi = wr*data[j+1]+wi*data[j];
                data[j] = data[i]-tempr;
                data[j+1] = data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }

            wtemp = wr;
            wr = wr*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }

        mmax = istep;
    }

}
