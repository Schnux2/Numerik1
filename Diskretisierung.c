#include <stdio.h>
#include <math.h>
#include "Bisektion_Newton.h"
#include "Diskretisierung.h"

/*
Diskretisierung des Randwertproblems u''(x)=f(u(x)) mit den Startwerten (x0, y0) und (x1,y1) (wobei x1>x0).
Für die Approximation wird auch noch f' benötigt.
Dafür wird das Newton-Verfahren und eine Gauß-Elimination angewandt. Der Vektor u hat die Länge (x1-x0)*ANZAHL_WERTE und
enthält die Funktionswerte an dieser Position. iterNum ist die Anzahl der Iterationsschritte
Das Ergebnis wird in der Datei "disk.txt" ausgegeben, in der Gestalt von (x,y)-Wertepaaren. Wird 0 zurückgegeben, so war die
Funktion erfolgreich, sonst wird 1 zurückgegeben.
*/
int Randwertproblem(double (*f)(double x), double (*df)(double x), double x0, double x1, double y0, double y1, double u[], int iterNum){
	if (x0>=x1){
		return 1;
	}
	double h=(x1-x0)/(ANZAHL_WERTE+1);
	//Startwerte, aber eher zufällige (geradlinige Verbindung)
	for (int i=0;i<ANZAHL_WERTE;i++){
		u[i]=(y1-y0)/(x1-x0)*i/(ANZAHL_WERTE+1);
	}
	
	//Es ist u^(k+1)=u^(k)+v^(k) eine bessere Approximation als u^(k) (vgl. Newton-Verfahren für Vektoren,
	//die zu betrachtende Matrix berechnet sich aus den zweiten Ableitungen).
	for (int k=0;k<iterNum;k++){
		double v[ANZAHL_WERTE];
		//Lösen der durch die Approximation der zweiten Ableitungen entstehenden linearen Gleichung für v mit einer LR-Zerlegung
		double r[ANZAHL_WERTE];
		//double l[ANZAHL_WERTE];
		r[0]=2+h*h*df(u[0]);
		//l[i]=-1/r[i-1]
		for (int i=1;i<ANZAHL_WERTE;i++){
			r[i]=2+h*h*df(u[i])-1/r[i-1];
		}
		
		//"L"-Schritt
		v[0]=y0-2*u[0]+u[1]-h*h*f(u[0]);
		for (int i=1;i<ANZAHL_WERTE-1;i++){
			v[i]=u[i-1]-2*u[i]+u[i+1]-h*h*f(u[i])+1/(r[i-1])*v[i-1];
		}
		v[ANZAHL_WERTE-1]=u[ANZAHL_WERTE-2]-2*u[ANZAHL_WERTE-1]+y1-h*h*df(u[ANZAHL_WERTE-1])+1/r[ANZAHL_WERTE-2]*v[ANZAHL_WERTE-2];
		
		//"R"-Schritt
		v[ANZAHL_WERTE-1]=1/r[ANZAHL_WERTE-1]*v[ANZAHL_WERTE-1];
		for (int i=ANZAHL_WERTE-2;i>-1;i--){
			v[i]=1/r[i]*(v[i]+v[i+1]);
		}
		
		//Addition von v zu u
		for (int i=0;i<ANZAHL_WERTE;i++){
			u[i]=u[i]+v[i];
		}
		
		//Abbruch, wenn genau genug
		double maxChange=0;
		for (int i=0;i<ANZAHL_WERTE;i++){
			if (fabs(v[i])>maxChange){
				maxChange=fabs(v[i]);
			}
		}
		if (maxChange<GENAUIGKEIT){
			break;
		}
	}
	
	
	FILE* dat=fopen("./disk.txt","w");
	for (int i=0;i<ANZAHL_WERTE;i++){
		fprintf(dat,"%f\t%f\n",x0+i*(x1-x0)/(ANZAHL_WERTE-1),u[i]);
	}
	fclose(dat);
	return 0;
}