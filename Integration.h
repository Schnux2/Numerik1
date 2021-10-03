#ifndef INTEGRATION_H_
#define INTEGRATION_H_

double Newton_Cotes(double (*f)(double x), double a, double b, int n);
double Newton_Cotes_Iteriert(double (*f)(double x), double a, double b, int n, int teil_Anzahl);
double Romberg(double (*f)(double x), double a, double b, int min_exp, int max_deg);

#endif