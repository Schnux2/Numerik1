#include <math.h>
#include "matrix.h"
#include "LR_Zerlegung.h"

/*
LR-Zerlegung mit Spaltenpivotsuche der Matrix a, falls möglich. Das Resultat wird wieder in a abgespeichert,
da bekannt ist, dass in der Diagonale der L-Matrix nur Einsen stehen. perm ist ein Zeiger, in dem die verwendete Permutation
abgespeichert wird.
Der Rückgabewert ist 0 bei Erfolg und 1 sonst.
*/
int LR(int dim, double* a[], int perm[]){
	for (int k=0;k<dim;k++){
		perm[k]=k;
	}
	for (int k=0;k<dim-1;k++){
		int p=k;
		for (int i=k+1;i<dim;i++){
			if (fabs(a[i][k])>fabs(a[p][k])){
				p=i;
			}
		}
		if (p!=k){
			int dum=perm[k];
			perm[k]=perm[p];
			perm[p]=dum;
			for (int j=0;j<dim;j++){
				double dum=a[p][j];
				a[p][j]=a[k][j];
				a[k][j]=dum;
			}
		}
		for (int i=k+1;i<dim;i++){
			if (a[k][k]!=0){
				a[i][k]=a[i][k]/a[k][k];
			}
			else{
				return 1; //LR-Zerlegung nicht möglich, da Matrix nicht invertierbar
			}
			for (int j=k+1;j<dim;j++){
				a[i][j]=a[i][j]-a[i][k]*a[k][j];
			}
		}
	}
	return 0;
}


/*
Vorwärtssubstitution, d.h. Lösung der Gleichung L*x=b, wobei L eine untere Dreiecksmatrix ist, bei der in der Diagonalen nur Einsen stehen.
Das Resultat wird in b abgespeichert. Der Rückgabewert ist 0 bei Erfolg (also eigentlich immer).
*/
int vorwaerts(int dim, double* L[], double b[]){
	for (int i=0;i<dim;i++){
		double h=0;
		for (int j=0;j<i;j++){
			h=h+L[i][j]*b[j];
		}
		b[i]=(b[i]-h)/1; // /L[i][i] im allgemeinen Fall
	}
	return 0;
}

/*
Rückwärtssubstitution, d.h. Lösung der Gleichung R*x=b, wobei R eine obere Dreiecksmatrix ist
Das Resultat wird in b abgespeichert. Der Rückgabewert ist 0 bei Erfolg und 1 sonst.
*/
int rueckwaerts(int dim, double* R[], double b[]){
	for (int i=dim-1;i>=0;i--){
		double h=0;
		for (int j=dim-1;j>i;j--){
			h=h+R[i][j]*b[j];
		}
		if (R[i][i]==0){
			return 1;
		}
		b[i]=(b[i]-h)/R[i][i];
	}
	return 0;
}

/*
Löst das Gleichungssystem A+x=b mit einer LR-Zerlegung. Der Rückgabwert ist 0 bei Erfolg und 1 sonst.
Das Ergebnis wird in b abgespeichert und auch A wird verändert.
*/
int loese(int dim, double* a[], double b[]){
	int perm[dim];
	int success=LR(dim,a,perm);
	
	if (success!=0){
		return 1;
	}
	
	//Aus Ax=b folgt PAx=Pb => Löse letztere Gleichung, für die die LR-Zerlegung vorhanden ist
	double c[dim];
	for (int i=0;i<dim;i++){
		c[i]=b[perm[i]];
	}
	
	success=vorwaerts(dim,a,c);
	if (success!=0){
		return 1;
	}
	success=rueckwaerts(dim,a,c);
	if (success!=0){
		return 1;
	}
	
	for (int i=0;i<dim;i++){
		b[i]=c[i];
	}
	return 0;
}