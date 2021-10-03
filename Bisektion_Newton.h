#ifndef BISEKTION_NEWTON_H_
#define BISEKTION_NEWTON_H_

#define GENAUIGKEIT 0.001
#define REK_MAX 100

double Bisektion(double (*f)(double x), double li, double re, int rekTiefe);
double Newton(double (*f)(double x), double (*df)(double x), double x, int rekTiefe);

#endif