#ifndef DISKRETISIERUNG_H_
#define DISKRETISIERUNG_H_

#define ANZAHL_WERTE 1024

int Randwertproblem(double (*d2u)(double x), double (*d3u)(double x), double x0, double x1, double y0, double y1, double u[], int iterNum);

#endif