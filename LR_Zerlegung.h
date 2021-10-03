#ifndef LR_ZERLEGUNG_H_
#define LR_ZERLEGUNG_H_

int LR(int dim, double* a[], int perm[]);
int vorwaerts(int dim, double* L[], double b[]);
int rueckwaerts(int dim, double* R[], double b[]);
int loese(int dim, double* a[], double b[]);
#endif