#include <math.h>
#include "Interpolation.h"

/*
Interpoliert die Punkte (xi,fi) mit einem Polynom vom Grad n (wobei n+1 die Länge von xi, fi und koeff sein muss
und die xi paarweise verschieden).
Die Koeffizienten werden in koeff abgespeichert. Der Rückgabewert ist 0 bei Erfolg und 1 sonst.
Es wird die Newton-Basis und zur Berechnung die Newtonsche Interpolationsformel verwendet (Dreiecksschema).
*/
int Interpolation(int n, double xi[], double fi[], double koeff[]){
	for (int i=0;i<n+1;i++){
		for (int j=0;j<n+1;j++){
			if (i!=j && xi[i]==xi[j]){
				return 1;
			}
		}
	}
	
	for (int i=0;i<n+1;i++){
		koeff[i]=fi[i];
	}
	
	//Berechnung "von links nach rechts", wobei immer nur die nötigen Einträge gespeichert werden
	for (int i=1;i<n+1;i++){
		//Berechnung der i-ten Spalte im Dreiecksschema "von unten nach oben", spart einen Speicherplatz
		for (int j=n;j>i-1;j--){
			koeff[j]=(koeff[j]-koeff[j-1])/(xi[j]-xi[j-i]);
		}
	}
	return 1;
}

/*
Berechnung eines gegebenen Interpolationspolynoms mit den Stützstellen xi und den Koeffizienten in der Newton-Basis
di mittels des Horner-Schemas (Rückwärtsrekursion) im Punkt x. w ist der im Horner-Schema zu übergebende Wert,
zu Beginn 0. k ist zu Beginn gleich n+1.
*/
double Horner(int k, double x, double w, double xi[], double di[]){
	if (k<=-1){
		return w;
	}
	return Horner(k-1,x,di[k]+w*(x-xi[k]),xi,di);
}

/*
Berechnung des Werts eines die Punkte (xi,fi) interpolierenden Polynoms im Punkt x mit dem Neville/Aitken-Schema. 
*/
double Neville(int n, double x, double xi[], double fi[]){
	for (int i=0;i<n+1;i++){
		for (int j=0;j<n+1;j++){
			if (i!=j && xi[i]==xi[j]){
				return NAN;
			}
		}
	}
	
	double werte[n+1];
	for (int i=0;i<n+1;i++){
		werte[i]=fi[i];
	}
	
	//Berechnung "von links nach rechts", wobei immer nur die nötigen Einträge gespeichert werden, ähnlich wie
	//bei der Berechnung des Interpolationspolynoms
	for (int i=1;i<n+1;i++){
		//Berechnung der i-ten Spalte im Dreiecksschema "von unten nach oben"
		for (int j=n;j>i-1;j--){
			werte[j]=werte[j]+(xi[j]-x)/(xi[j]-xi[j-i])*(werte[j-1]-werte[j]);
		}
	}
	return werte[n];
}