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

    // S_1
    double *S1 = malloc(N*N*sizeof(double));

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

    // compute S_1
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            // wat
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



