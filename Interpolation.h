#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

int Interpolation(int n, double xi[], double fi[], double koeff[]);
double Horner(int k, double x, double w, double xi[], double ai[]);
double Neville(int n, double x, double xi[], double fi[]);

#endif