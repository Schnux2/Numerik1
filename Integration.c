#include <math.h>
#include <stdio.h>
#include "Integration.h"
#include "matrix.h"
#include "LR_Zerlegung.h"

/*
Numerische Integration der gegebenen Funktion mit der Newton-Cotes-Formel der Ordnung n von a bis b
*/
double Newton_Cotes(double (*f)(double x), double a, double b, int n){
	double xi[n+1]; //äquidistante Stützstellen in der Newton-Cotes-Formel
	for (int i=0;i<n+1;i++){
		xi[i]=a+(b-a)*i/n;
	}
	
	//Berechnung der Koeffizienten der Koeffizienten für die Quadraturformel
	double** mat=allocM(n+1,n+1);
	double koeff[n+1];
	for (int i=0;i<n+1;i++){
		for (int j=0;j<n+1;j++){
			mat[i][j]=pow(xi[j],i);
		}
		koeff[i]=(pow(b,i+1)-pow(a,i+1))/(i+1);
	}
	loese(n+1,mat,koeff);
	freeM(n+1,n+1,mat);
	
	//Approximation des Integrals mit der Newton-Cotes-Formel
	double res=0;
	for (int i=0;i<n+1;i++){
		res=res+koeff[i]*f(xi[i]);
	}
	return res;
}

/*
Numerische Berechnung eines Integrals mithilfe einer iterierten Newton-Cotes-Formel, wobei das Integrationsintervall
in teil_Anzahl Teile geteilt wird.
*/
double Newton_Cotes_Iteriert(double (*f)(double x), double a, double b, int n, int teil_Anzahl){
	double res=0;
	for (int k=0;k<teil_Anzahl;k++){
		double ak=a+k*(b-a)/teil_Anzahl;
		double bk=a+(k+1)*(b-a)/teil_Anzahl;
		res=res+Newton_Cotes(f,ak,bk,n);
	}
	return res;
}

/*
Numerische Berechnung eines Integrals mithilfe einer Romberg-Folge. Hierbei wird zunächst festgestellt, dass eine
iterierte Trapezregel mit kleiner werdenden Teilintervallbreiten h das Integral immer genauer approximieren kann
(Euler-MacLaurinsche Summenformel).
Auf diese Weise enthält man eine Funktion in Abhängigkeit von h, deren Wert bei 0 gleich dem Integral sein muss.
Diese Funktion kann man nun mit einem Neville-Aitken-Schema interpolieren. Hierbei werden für die Romberg-Folge
zunächst die Approximationen für die iterierte Trapezregel mit h=2^k*(a-b) berechnet (wobei k<0 ist). Der minimale
Exponent ist hier mit min_exp angegeben.
Allerdings wird für die Euler-MacLaurinsche Summenformel eine hohe Differenzierbarkeitsordnung benötigt, die häufig nicht
gegeben ist. Aus diesem Grund kann noch die maximale Ordnung des Interpolationspolynoms max_deg angegeben werden,
wodurch die Approximation besser möglich ist (das Dreiecks-Schema von Neville-Aitken (das Extrapolationsschema))
wird dann nicht bis "ganz rechts" durchgeführt).
*/
double Romberg(double (*f)(double x), double a, double b, int min_exp, int max_deg){
	if (min_exp>0){
		return NAN;
	}
	
	//Berechnung des Dreiecksschemas "von oben nach unten" (so kann man leichter ein neues
	//Romberg-Folgenglied dranhängen)
	//Gespeichert wird nur die unterste Zeile des Dreiecksschemas.
	double werte[-min_exp<max_deg ? -min_exp : max_deg];
	for (int k=0;k>min_exp-1;k--){
		double dum=werte[0];
		werte[0]=Newton_Cotes_Iteriert(f,a,b,1,pow(2,-k));
		printf("%f\t",werte[0]);
		for (int i=1;i<(-k<max_deg ? -k : max_deg)+1;i++){
			double dum2=werte[i];
			werte[i]=werte[i-1]+(werte[i-1]-dum)/(pow(4,i)-1);
			dum=dum2;
			printf("%f\t",werte[i]);
		}
		printf("\n");
	}
	return werte[(-min_exp<max_deg ? -min_exp : max_deg)-1];
}