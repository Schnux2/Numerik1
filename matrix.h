#ifndef MATRIX_H_
#define MATRIX_H_
//wird nur einmal eingebunden

void matProd(int m, int n, int p, double** a, double** b, double** ret);
void matVecProd(int m, int n, double** a, double b[], double ret[]);
double det(int m, double** a);
double** copyM(int m, int n, double** a);
double** allocM(int m, int n);
void freeM(int m, int n, double* mat[]);
void printM(int m, int n, double* mat[]);
void printV(int n, double v[]);
void vecAdd(int n, double a[], double b[], double res[]);
double vecLen(int n, double v[]);

#endif