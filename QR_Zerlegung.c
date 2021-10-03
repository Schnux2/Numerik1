#include <math.h>
#include "matrix.h"
#include "LR_Zerlegung.h"
#include "QR_Zerlegung.h"

/*
Berechnung der QR-Zerlegung der Matrix a mit Householder-Matrizen (Hyperebenenspiegelungen). Die R-Matrix wird in a abgespeichert,
sowie die Vektoren der Hyperebenenspiegelung unterhalb der Diagonale. Die Diagonaleinträge von R werden
gesondert in diag abgespeichert. Der Rückgabewert ist 0 bei Erfolg und 1 sonst (also eigentlich immer 0).
*/
int QR(int m, int n, double* a[], double diag[]){
	for (int k=0;k< (n>m ? m : n);k++){
		diag[k]=a[k][k]*a[k][k];
		for (int j=k+1;j<m;j++){
			diag[k]=diag[k]+a[j][k]*a[j][k];
		}
		diag[k]=-(a[k][k]!=0 ? a[k][k]/fabs(a[k][k]) : 1) *sqrt(diag[k]);
		
		double b=diag[k]*(a[k][k]-diag[k]);
		a[k][k]=a[k][k]-diag[k];
		for (int i=k+1;i<n;i++){
			double c=a[k][k]*a[k][i];
			for (int j=k+1;j<m;j++){
				c=c+a[j][k]*a[j][i];
			}
			double d=c/b;
			for (int j=k;j<m;j++){
				a[j][i]=a[j][i]+d*a[j][k];
			}
		}
	}
	return 0;
}

/*
Minimiert |a*x-b|^2 mithilfe einer QR-Zerlegung und gibt das Ergebnis in x zurück (lineares Ausgleichsproblem).
Die Dimension der Matrix ist m*n. Das Resultat wird wieder in b zurückgegeben (auch wenn b eigentlich zu lang ist).
Ist die Methode nicht anwendbar (z.B. mit n>m), so wird 1 zurückgegeben, sonst 0.
Der Algorithmus berechnet zunächst die QR-Zerlegung von a und dann Q^t*b. Dann nimmt R die Gestalt
(R1) an und der zu minimierende Betrag berechnet sich als |(R1)*x - (b1)|, wobei Q^t*b=(b1).
(0 )                                                      |(0 )     (b2)|              (b2)
Das ist der Fall, wenn R1*x=b1 bzw. x=R1^{-1}*b1. Dieses Gleichungssystem in Dreiecksgestalt kann einfach gelöst werden
(wobei davon ausgegangen wird, dass R1 vollen Rang hat).
*/
int Ausgleichsproblem(int m, int n, double* a[], double b[]){
	if (n>m){
		return 1;
	}
	double diag[n];
	QR(m,n,a,diag);
	
	//Berechnung von Q^t*b. Es ist Q=Q^(n)*...*Q^(1), wobei Q^(k) orthogonale Householder-Matrizen sind, d.h. von der
	//Form I_m-2*v^t*v.
	
	//Durchlaufen der Householder-Matrizen
	for (int k=0;k<n;k++){
		double c[m]; //Ergebnis der Multiplikation, wird dann wieder in b gespeichert
		
		double v[m]; //Spiegelungsvektor der Householder-Matrix
		for (int i=0;i<k;i++){
			v[i]=0;
		}
		for (int i=k;i<m;i++){
			v[i]=a[i][k];
		}
		
		//Reskalieren des Spiegelungsvektors
		double len=vecLen(m,v);
		for (int j=0;j<m;j++){
			v[j]=v[j]/len;
		}
		
		//Bilden der Matrix-Vektor-Produkts, wobei die Struktur der Householder-Matrizen berücksichtigt wird
		for (int i=0;i<m;i++){
			c[i]=b[i];
			for (int j=0;j<m;j++){
				c[i]=c[i]-2*v[i]*v[j]*b[j];
			}
		}
		
		//Kopieren von c in b, für den nächsten Iterationsschritt
		for (int i=0;i<m;i++){
			b[i]=c[i];
		}
	}
	
	//damit rueckwaerts angewandt werden kann
	for (int j=0;j<n;j++){
		a[j][j]=diag[j];
	}
	
	int success=rueckwaerts(n,a,b); //b ist länger als nötig und a hat z.T. falsche Einträge, aber das ist egal
	if (success!=0){
		return 1;
	}
	
	return 0;
}