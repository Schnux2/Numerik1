#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "LR_Zerlegung.h"
#include "QR_Zerlegung.h"
#include "Bisektion_Newton.h"
#include "Diskretisierung.h"
#include "Interpolation.h"
#include "Integration.h"
#include "SOR.h"

int testMat();
int testLR();
int testQR();
int testQR2();
int testBisektion();
int testDisk();
int testInter();
int testInt();
int testSOR();

double f(double x);
double df(double x);
double g(double x);
double h(double x, double y);
double k(double x, double y);

int main(int argc, char* argv[]){
	testSOR();
	return 0;
}

int testSOR(){
	for (int n=10;n<50;n=n+10){
		double exact[(n-1)*(n-1)]; //exakte Lösung für f=-4, g=x^2+y^2
		for (int i=0;i<(n-1)*(n-1);i++){
			exact[i]=((double) ((i/(n-1)+1)*(i/(n-1)+1)+(i%(n-1)+1)*(i%(n-1)+1)))/(n*n);
		}
		
		for (int omega10=5;omega10<20;omega10++){
			//Startwerte
			double u[(n-1)*(n-1)];
			for (int i=0;i<(n-1)*(n-1);i++){
				u[i]=0;
			}
			double omega=((double) omega10)/10;
			int iterNum=0;
			double diffToExact=0;
			for (int i=0;i<(n-1)*(n-1);i++){
				diffToExact=diffToExact+(u[i]-exact[i])*(u[i]-exact[i]);
			}
			diffToExact=sqrt(diffToExact);
			double diffToExactStart=diffToExact;
			
			//Iteration so lange, bis die Differenz zum exakten Wert kleiner als 0.1*der Startdifferenz ist
			while (diffToExact>0.1*diffToExactStart && iterNum<50){
				iterNum++;
				Poissonproblem_2D(h,k,n,1,u,omega);
				diffToExact=0;
				for (int i=0;i<(n-1)*(n-1);i++){
					diffToExact=diffToExact+(u[i]-exact[i])*(u[i]-exact[i]);
				}
				
				diffToExact=sqrt(diffToExact);
			}
			printf("Das SOR-Verfahren mit n=%d, omega=%f benötigte %d Iterationsschritte.\n", n, omega, iterNum);
		}
		
		//optimaler Relaxationsparameter omega (ohne Begründung)
		double u[(n-1)*(n-1)];
		for (int i=0;i<(n-1)*(n-1);i++){
			u[i]=0;
		}
		double omega=2/(1+sin(acos(-1)/n));
		int iterNum=0;
		double diffToExact=0;
		for (int i=0;i<(n-1)*(n-1);i++){
			diffToExact=diffToExact+(u[i]-exact[i])*(u[i]-exact[i]);
		}
		diffToExact=sqrt(diffToExact);
		double diffToExactStart=diffToExact;
		
		while (diffToExact>0.1*diffToExactStart && iterNum<50){
			iterNum++;
			Poissonproblem_2D(h,k,n,1,u,omega);
			diffToExact=0;
			for (int i=0;i<(n-1)*(n-1);i++){
				diffToExact=diffToExact+(u[i]-exact[i])*(u[i]-exact[i]);
			}
			
			diffToExact=sqrt(diffToExact);
		}
		printf("Das SOR-Verfahren mit n=%d, omega=%f benötigte %d Iterationsschritte.\n", n, omega, iterNum);
	}
	
	/*
	int n=8;
	double u[(n-1)*(n-1)];
	for (int i=0;i<(n-1)*(n-1);i++){
		u[i]=0;
	}
	double omega=1;
	int iterNum=10;
	Poissonproblem_2D(h,k,n,iterNum,u,omega);*/
	return 0;
}

int testInt(){
	double r=Romberg(g,0,1,-4,3);
	printf("Numerische Integration: %f",r);
	return 0;
}

int testInter(){
	int n=5;
	double xi[n+1];
	double fi[n+1];
	double koeff[n+1];
	for (int i=0;i<n+1;i++){
		xi[i]=-5+((double) i)*10/n;
		fi[i]=1/(1+xi[i]*xi[i]);
	}
	
	Interpolation(n,xi,fi,koeff);
	printf("Interpolationspolynom:\n");
	for (int i=0;i<n+1;i++){
		printf("%f",koeff[i]);
		for (int j=0;j<i;j++){
			printf("*(x-%f)",xi[j]);
		}
		if (i<n){
			printf(" + ");
		}
	}
	printf("\n");
	
	//Tschebyscheff-Stützstellen
	for (int i=0;i<n+1;i++){
		xi[i]=0+(5+5)/2*cos((2*((double)i)+1)*3.14/(2*n+2));
		fi[i]=1/(1+xi[i]*xi[i]);
	}
	
	Interpolation(n,xi,fi,koeff);
	printf("Interpolationspolynom bzgl. Tschebyscheff-Stützstellen:\n");
	for (int i=0;i<n+1;i++){
		printf("%f",koeff[i]);
		for (int j=0;j<i;j++){
			printf("*(x-%f)",xi[j]);
		}
		if (i<n){
			printf(" + ");
		}
	}
	printf("\n");
	
	clock_t tstart, tend;
	int m=100;
	double xi2[m+1];
	double fi2[m+1];
	for (int i=0;i<m+1;i++){
		xi2[i]=-5+((double) i)*10/m;
		fi2[i]=1/(1+xi[i]*xi[i]);
	}
	tstart=clock();
	for (int k=0;k<10000;k++){
		Neville(m,0,xi2,fi2);
	}
	tend=clock();
	printf("Das Neville-Schema benötigte %f Sekunden.\n",((double) (tend-tstart))/CLOCKS_PER_SEC);
	
	tstart=clock();
	double koeff2[m+1];
	Interpolation(m,xi2,fi2,koeff2);
	for (int k=0;k<10000;k++){
		Horner(m+1,0,0,xi2,koeff2);
	}
	tend=clock();
	printf("Das Horner-Schema benötigte %f Sekunden.\n",((double) (tend-tstart))/CLOCKS_PER_SEC);
	return 0;
}

int testDisk(){
	double x0=0;
	double x1=1;
	double y0=4;
	double y1=1;
	double u[ANZAHL_WERTE];
	Randwertproblem(f, df, x0, x1, y0, y1, u, 10);
	return 0;
}

int testBisektion(){
	double n=Newton(sin,cos,2,0);
	printf("%f",n);
	return 0;
}

int testQR(){
	int m=3;
	int n=3;
	double** a=allocM(m,n);
	double b[m];
	
	a[0][0]=3; a[0][1]=1; a[0][2]=2;
	a[1][0]=1; a[1][1]=5; a[1][2]=1;
	a[2][0]=2; a[2][1]=1; a[2][2]=4;
	
	b[0]=3; b[1]=-4; b[2]=-1;
	
	Ausgleichsproblem(m,n,a,b);
	printV(n,b);
	
	freeM(m,n,a);
	
	/*
	int m=5;
	int n=3;
	double** a=allocM(m,n);
	double diag[n];
	
	a[0][0]=1; a[0][1]=-2; a[0][2]=1;
	a[1][0]=1; a[1][1]=3; a[1][2]=-1;
	a[2][0]=1; a[2][1]=0; a[2][2]=-1;
	a[3][0]=1; a[3][1]=3; a[3][2]=-1;
	a[4][0]=0; a[4][1]=sqrt(7); a[4][2]=11/sqrt(7);
	
	QR(m,n,a,diag);
	printM(m,n,a);
	printV(n,diag);
	
	freeM(m,n,a);*/
	return 0;
}

int testQR2(){
	double Staerke[] = {6,75,87,55,34,98,91,45,51,17,36,97,74,24,85,96,92,94,84,99};
	double Protein[] = {10.3,12.2,14.5,11.1,10.9,18.1,14.0,10.8,11.4,11.0,10.2,17.0,13.8,10.1,14.4,15.8,15.6,15.0,13.3,19.0};
	double** werte=allocM(20,3);
	for (int i=0;i<20;i++){
		werte[i][0]=Staerke[i]*Staerke[i];
		werte[i][1]=Staerke[i];
		werte[i][2]=1;
	}
	printM(20,3,werte);
	
	Ausgleichsproblem(20,3,werte,Protein);
	freeM(20,3,werte);
	printf("Die am besten geeignete Funktion ist %f*x^2+%f*x+%f.\n",Protein[0],Protein[1],Protein[2]);
	return 0;
}

int testLR(){
	int dim=3;
	double b[dim];
	double** mat=allocM(dim,dim);
	
	mat[0][0]=2; mat[0][1]=1; mat[0][2]=1;
	mat[1][0]=4; mat[1][1]=2; mat[1][2]=-1;
	mat[2][0]=-1; mat[2][1]=0; mat[2][2]=7;
	
	b[0]=1; b[1]=2; b[2]=3;
	
	printf("Löse \n(");
	for (int i=0;i<dim;i++){
		for (int j=0;j<dim;j++){
			printf("%f", mat[i][j]);
			if (j<dim-1){
				printf(", ");
			}
			else{
				printf(")");
				if (i==0){
					printf(" * x =\t(%f)\n",b[i]);
				}
				else{
					printf("\t\t(%f)\n",b[i]);
				}
			}
		}
		if (i<dim-1){
			printf("(");
		}
	}
	printf("\n");
	
	loese(dim,mat,b);
	
	for (int i=0;i<dim;i++){
		if (i==0){
			printf("x = ");
		}
		printf("\t(%f)\n",b[i]);
	}
	
	freeM(dim,dim,mat);
	return 0;
}

int testMat(){
	double b[3];
	double** mat=allocM(2,3);
	double** mat2=allocM(3,4);
	
	mat[0][0]=2; mat[0][1]=1; mat[0][2]=1;
	mat[1][0]=4; mat[1][1]=2; mat[1][2]=-1;
	
	mat2[0][0]=2; mat2[0][1]=1; mat2[0][2]=1; mat2[0][3]=4;
	mat2[1][0]=4; mat2[1][1]=2; mat2[1][2]=-1; mat2[1][3]=2;
	mat2[2][0]=2; mat2[2][1]=1; mat2[2][2]=1; mat2[2][3]=-2;
	
	b[0]=1; b[1]=2; b[2]=3;
	
	printM(2,3,mat);
	printM(3,4,mat2);
	printV(3,b);
	
	double** mat3=allocM(2,4);
	matProd(2,3,4,mat,mat2,mat3);
	printM(2,4,mat3);
	double c[2];
	matVecProd(2,3,mat,b,c);
	printV(2,c);
	printf("%f\n",vecLen(2,c));
	
	double** mat4=allocM(3,3);
	mat4[0][0]=2; mat4[0][1]=1; mat4[0][2]=1;
	mat4[1][0]=4; mat4[1][1]=2; mat4[1][2]=-1;
	mat4[2][0]=2; mat4[2][1]=2; mat4[2][2]=3;
	double** mat5=copyM(3,3,mat4);
	printM(3,3,mat5);
	double d=det(3,mat5);
	printf("%f\n",d);
	freeM(2,3,mat);
	freeM(3,4,mat2);
	freeM(2,4,mat3);
	freeM(3,3,mat4);
	freeM(3,3,mat5);
	return 0;
}


double f(double x){
	return x*x*5/2;
}
double df(double x){
	return 5*x;
}
double g(double x){
	return (16*x-16)/(pow(x,4)-2*pow(x,3)+4*x-4);
}
double h(double x, double y){
	return -4;
}
double k(double x, double y){
	return x*x+y*y;
}