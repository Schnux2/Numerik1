#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

void matProd(int m, int n, int p, double** a, double** b, double** ret){
	for (int i=0;i<m;i++){
		for (int j=0;j<p;j++){
			ret[i][j]=0;
			for (int k=0;k<n;k++){
				ret[i][j]=ret[i][j]+a[i][k]*b[k][j];
			}
		}
	}
}

void matVecProd(int m, int n, double** a, double b[], double ret[]){
	for (int i=0;i<m;i++){
		ret[i]=0;
		for (int j=0;j<n;j++){
			ret[i]=ret[i]+a[i][j]*b[j];
		}
	}
}

double det(int m, double** a){
	if (m<=1){
		return a[0][0];
	}
	double res=0;
	for (int j=0;j<m;j++){
		double** sub=allocM(m-1,m-1);
		for (int i=0;i<m-1;i++){
			for (int k=0;k<m;k++){
				int k2=k;
				if (k>=j){
					k2=k2+1;
				}
				sub[i][k]=a[i+1][k2];
			}
		}
		res=res+pow(-1,j)*a[0][j]*det(m-1,sub);
		free(sub);
	}
	return res;
}

double** copyM(int m, int n, double** a){
	double** b=allocM(m,n);
	for (int i=0;i<m;i++){
		for (int j=0;j<n;j++){
			b[i][j]=a[i][j];
		}
	}
	return b;
}

double** allocM(int m, int n){
	double** mat=malloc(sizeof(double)*m*n);
	for (int i=0;i<m;i++){
		mat[i]=malloc(sizeof(double)*n);
	}
	return mat;
}

void freeM(int m, int n, double* mat[]){
	for (int i=0;i<m;i++){
		free(mat[i]);
	}
	free(mat);
}

void printM(int m, int n, double* mat[]){
	for (int i=0;i<m;i++){
		printf("(");
		for (int j=0;j<n;j++){
			printf("%f",mat[i][j]);
			if (j<n-1){
				printf(",\t");
			}
		}
		printf(")\n");
	}
}

void printV(int n, double v[]){
	for (int i=0;i<n;i++){
		printf("(%f)\n",v[i]);
	}
}

void vecAdd(int n, double a[], double b[], double res[]){
	for (int i=0;i<n;i++){
		res[i]=a[i]+b[i];
	}
}

double vecLen(int n, double v[]){
	double res=0;
	for (int i=0;i<n;i++){
		res=res+v[i]*v[i];
	}
	return sqrt(res);
}