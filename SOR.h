#ifndef SOR_H_
#define SOR_H_

int Poissonproblem_2D(double (*f)(double x, double y), double (*g)(double x, double y), int n, int iterNum, double u[], double omega);

#endif