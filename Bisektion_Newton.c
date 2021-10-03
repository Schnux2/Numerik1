#include <math.h>
#include "Bisektion_Newton.h"

/*
Bisektions-Verfahren zur Nullstellenberechnung, rekursiv. Wird keine Nullstelle gefunden, wird NaN zurückgegeben.
li und re sind die rechte und linke Grenze, rekTiefe ist die Tiefe der Rekursion (für eine Abbruchbedingung).
*/
double Bisektion(double (*f)(double x), double li, double re, int rekTiefe){
	double mi=(re+li)/2;
	if (fabs(f(mi))<GENAUIGKEIT){
		return mi;
	}
	if (rekTiefe>REK_MAX){
		return NAN;
	}
	rekTiefe++;
	if ((f(mi)>0 && f(li)<0) || (f(mi)<0 && f(li)>0)){
		return Bisektion(f,li,mi,rekTiefe);
	}
	else{
		return Bisektion(f,mi,re,rekTiefe);
	}
	return NAN; //wird nie erreicht
}

/*
Newton-Verfahren zur Nullstellenberechnung, rekursiv. Wird keine Nullstelle gefunden, wird NaN zurückgegeben.
f ist die Funktion, df die Ableitung, x der Startwert, rekTiefe und die Tiefe der Rekursion.
*/
double Newton(double (*f)(double x), double (*df)(double x), double x, int rekTiefe){
	if (fabs(f(x))<GENAUIGKEIT){
		return x;
	}
	if (rekTiefe>REK_MAX){
		return NAN;
	}
	rekTiefe++;
	double x1=x-f(x)/df(x);
	return Newton(f,df,x1,rekTiefe);
}