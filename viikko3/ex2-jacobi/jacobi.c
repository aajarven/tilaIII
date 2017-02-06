#include <math.h>


struct matrixIndex{
    int i;
    int j;
};


/**
 * Jacobi method as presented in http://www.cmi.ac.in/~ksutar/NLA2013/iterativemethods.pdf
**/

double* jacobi(double *Q, int N){
    // copy of Q
    double *D = malloc(N*N*sizeof(double));
    memcpy(D, Q, N*N*sizeof(double));

    // unit matrix
    double *S = malloc(N*N*sizeof(double));
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            if (i==j){
                S[i][j] = 1;
            } else {
                S[i][j] = 0;
            }
        }
    }

    // S_1 & S_1t
    double *S1 = malloc(N*N*sizeof(double));
    double *S1t = malloc(N*N*sizeof(double));

    // find biggest off diagonal value
    struct matrixIndex biggest = biggestOffDiag(D);

    // finding the rotational angle
    double theta;
    if (D[biggest.i][biggest.i] == D[biggest.j][biggest.j]){
        if (D[biggest.i][biggest.j] > 0){
            theta = PI/2;
        } else {
            theta = -PI/2;
        }
    } else {
        theta = 0.5*atan(2*D[biggest.i][biggest.j]/D[biggest.i][biggest.i]-D[biggest.j][biggest.j]);
    }

    // compute S_1 (and its transpose)
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            if (i==j){
                S1[i][j] = 1;
                S1t[i][j] = 1;
            } else {
                S1[i][j] = 0;
                S1t[i][j] = 0;
            }
        }
    }
    S1[biggest.i][biggest.i] = cos(theta);
    S1[biggest.j][biggest.j] = S1[biggest.i][biggest.i];
    S1[biggest.j][biggest.i] = sin(theta);
    S1[biggest.i][biggset.j] = -S1[biggest.j][biggest.i];
    S1t[biggest.i][biggest.i]=S1[biggest.i][biggest.i];
    S1t[biggest.j][biggest.j]=S1[biggest.j][biggest.j];
    S1t[biggest.i][biggest.j]=S1[biggest.j][biggest.i];
    S1t[biggest.j][biggest.i]=S1[biggest.i][biggest.j];

    // matrix for doing multiplication S1t*D*S1
    double *tmp = malloc(N*N*sizeof(double));

    // S1t*D
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            tmp[i][j] = 0;
            for (int k=0; k<N; k++){
                tmp[i][j] += S1t[i][k]*D[k][j];
            }
        }
    }

    // D*S1
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            D[i][j] = 0;
            for (int k=0; k<N; k++){
                D[i][j] += tmp[i][k]*S1[k][j];
            }
        }
    }

}

struct matrixIndex biggestOffDiag(double *matrix, int size){

    struct matrixIndex = {0, 0};

    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            if (fabs(matrix[i][j]) > fabs(matrix[matrixIndex.i][matrixIndex.j])){
                matrixIndex.i = i;
                matrixIndex.j = j;
            }
        }
    }

    return matrixIndex;



