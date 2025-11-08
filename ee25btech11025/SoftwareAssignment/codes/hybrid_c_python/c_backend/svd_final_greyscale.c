#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define epsilon 1e-12
#define iterations 20                       
//the no.of iterations that the code performs to efind each eigenvalue

//This function is used to normalize a vector preventing it from growing infinitely large
void normalize_vec(double *vec, int n) {
    double mag = 0;
    for (int i = 0; i < n; i++) mag += vec[i] * vec[i];
    mag = sqrt(mag);
    if (mag < epsilon) {                    //used to avoid dividing the vector by zero and growing it indefinitely large
        for (int i = 0; i < n; i++) vec[i] = 0.0;
        return;                 
    };
    for (int i = 0; i < n; i++) vec[i] /= mag;     //normalizing the vector
}

//This function is used to multiply a amtrix with a vector
void mult_mat_vec(double *mat, double *vec, int r, int c, double *res) {
    for (int i = 0; i < r; i++) {
        double sum = 0;
        for (int j = 0; j < c; j++)
            sum += mat[i * c + j] * vec[j];
        res[i] = sum;
    }
}

//The function uses iterative approach to find the largest eigenvalue and corresponsing eigenvector
double power_met(double *mat, int n, double *eigenvec) {
    for (int i = 0; i < n; i++)                             //creating a random with elements between 0 and 1
        eigenvec[i] = (double)rand() / RAND_MAX;

    double *temp = calloc(n, sizeof(double));               //creating a temporary array
    normalize_vec(eigenvec, n);

    for (int it = 0; it < iterations; it++) {               //iterative step to find the eigenvec
        mult_mat_vec(mat, eigenvec, n, n, temp);
        normalize_vec(temp, n);
        for (int j = 0; j < n; j++)
            eigenvec[j] = temp[j];
    }

    double lambda = 0;
    mult_mat_vec(mat, eigenvec, n, n, temp);
    for (int i = 0; i < n; i++)
        lambda += eigenvec[i] * temp[i];                    //loop to find the eigenvalue
    if (lambda < 0) lambda = 0;

    free(temp);
    return lambda;
}

//This function is used to find A transpose into A
void A_T_A(double *A, double *B, int r, int c) {
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++) {
            double sum = 0;
            for (int k = 0; k < r; k++)
                sum += A[k * c + i] * A[k * c + j];
            B[i * c + j] = sum;
        }
    }
}

//This function is used to capoy one matrix into another
void copy(double *B, double *B2, int n) {
    for (int i = 0; i < n; i++) B2[i] = B[i];
}

void find_svd(double *A, int r, int c, int k, double *f) {
    double *B = calloc(c * c, sizeof(double));
    double *B_def = calloc(c * c, sizeof(double));
    double *U = calloc(r * k, sizeof(double));
    double *V = calloc(c * k, sizeof(double));
    double *S = calloc(k, sizeof(double));
    double *eigenvec = calloc(c, sizeof(double));
    double *left_s_vec = calloc(r, sizeof(double));

    A_T_A(A, B, r, c);                                          //creating a symmetrix matrix B
    copy(B, B_def, c * c);     

    for (int i = 0; i < k; i++) {                               //loop used to find the first k singular and eigenvalues
        double lambda = power_met(B_def, c, eigenvec);
        S[i] = sqrt(lambda);
        for (int j = 0; j < c; j++)
            V[j * k + i] = eigenvec[j];                         //Storing right singular vectors in right singular matrix
        for (int p = 0; p < c; p++)                             //Deflating the matrix B to find the next largest eigenvalue
            for (int j = 0; j < c; j++)
                B_def[p * c + j] -= lambda * eigenvec[p] * eigenvec[j];
    }

    for (int i = 0; i < k; i++) {                               //Finding the left singular vectors and storing them in left singular matrix
        for (int j = 0; j < c; j++)
            eigenvec[j] = V[j * k + i];
        mult_mat_vec(A, eigenvec, r, c, left_s_vec);
        normalize_vec(left_s_vec, r); 
        for (int j = 0; j < r; j++)
            U[j * k + i] = left_s_vec[j];
    }

    for (int i = 0; i < r; i++) {                               //Reconstructing the final matrix
        for (int j = 0; j < c; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++)
                sum += U[i * k + p] * S[p] * V[j * k + p];
            f[i * c + j] = sum;
        }
    }

    free(B); free(B_def); free(U); free(V); free(S); free(eigenvec); free(left_s_vec);
}

//This function is used to frobenius norm.
double Frobenius(double *A, double *A_k, int r, int c){
    double e=0;
    for(int i = 0; i<r; i++){
        for(int j = 0; j<c; j++){
            double d = A[i*c +j]-A_k[i*c+j];
            e+=d*d;
        }
    }
    e= sqrt(e);
    return e;
}
