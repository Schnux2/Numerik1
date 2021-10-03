#include <stdio.h>
#include "SOR.h"

/*
Lösung des diskretisierten Poissonproblems auf (0,1)^2, -Laplace(u)=f auf (0,1)^2 und f=g auf dem Rand von (0,1)^2.
Dafür werden (n+1)^2 Gitterpunkte in diesen Bereich gelegt. Es ergibt sich ein zu lösendes (n-1)^2-dimensionales
Gleichungssystem (für die inneren Punkte). Dieses wird mit einem SOR-Verfahren mit dem Parameter omega genähert, wobei
iterNum Schritte durchgeführt werden.
Das Verfahren konvergiert, da die zu betrachtende Matrix das schwache Zeilensummenkriterium erfüllt und nicht zerfällt.
Die Gitterpunkte werden in u zeilenweise hintereinander abgespeichert (von unten nach oben, von links nach rechts).
Somit hat u die Länge (n-1)^2.
Die zu betrachtende Matrix ist von der Form A1=1/(h^2)*(T,  -Id, 0  ...   0) = L + D + R (L untere, R obere Dreiecksmatrix, D diagonal,
                                                       (-Id, T, -Id ...   0)
                                                       (0,  -Id, T  ...   0)
                                                       (0        ...      0)
                                                       (0     ...    -Id, T)
T=(4, -1, 0  ... 0)) und b=(b1    ) mit b1=(f(h,h)+(g(h,0)+g(0,h))/(h^2)    ), bj=(f(h,jh)+g(0,jh)/(h^2)  ) (für 1<j<n-1),
  (-1, 4, -1 ... 0)        (b2    )        (f(2h,h)+g(2h,0)/(h^2)           )     (f(2h,jh)               )
  (0, -1, 4  ... 0)        (...   )        (...                             )     (...                    )
  (0     ...    -1)        (...   )        (f(1-2h,h)+g(1-2h,0)/(h^2)       )     (f(1-2h,jh)             )
  (0  ...    -1, 4)        (b(n-1))        (f(1-h,h)+(g(1-h,0)+g(1,h))/(h^2))     (f(1-h,jh)+g(1,jh)/(h^2))
b(n-1)=(f(h,1-h)+(g(h,1)+g(0,1-h))/(h^2)    ) (wobei h=1/(n+1) und die zu lösende Gleichung A1*x=b ist).
       (f(2h,1-h)+g(2h,1)/(h^2)             )
       (...                                 )
       (f(1-2h,1-h)+g(1-2h,1)/(h^2)         )
       (f(1-h,1-h)+(g(1-h,1)+g(1,1-h))/(h^2))
Für das SOR-Verfahren wird eine Iteration x^(k+1) = (Id-B^(-1)*A1)*x^(k) + B^(-1)*b durchgeführt, wobei B=(1/omega*D+L).
Das bedeutet auch x^(k+1) = x^(k)-omega*D^(-1)*(L*x^(k+1)+(D+R)*x^(k)-b) oder für die einzelnen Komponenten
x^(k+1)_i = x^(k)_i - omega*1/a_ii*(\sum_{j=1}^{i-1} a_ij*x^(k+1)_j + \sum_{j=i}^n a_ij*x^(k)_j - b_i).
Das Ergebnis wird in der Datei poisson.txt abgespeichert. Der Rückgabewert ist immer 0.
*/
int Poissonproblem_2D(double (*f)(double x, double y), double (*g)(double x, double y), int n, int iterNum, double u[], double omega){
	double h=1/((double) n);
	double h2=h*h;
	
	//Anlegen von b
	//A1 wird aufgrund der einfachen Gestalt nicht als Matrix angelegt (bräuchte viel Speicherplatz)
	double b[(n-1)*(n-1)];
	b[0]=f(h,h)+(g(h,0)+g(0,h))/h2;
	for (int i=1;i<n-2;i++){
		b[i]=f(i*h,h)+g(i*h,0)/h2;
	}
	b[n-2]=f(1-h,h)+(g(1-h,0)+g(1,h))/h2;
	
	for (int j=1;j<n-2;j++){
		b[j*(n-1)]=f(h,j*h)+g(0,j*h)/h2;
		for (int i=1;i<n-2;i++){
			b[j*(n-1)+i]=f(i*h,j*h);
		}
		b[j*(n-1)+n-2]=f(1-h,j*h)+g(1,j*h)/h2;
	}
	
	b[(n-2)*(n-1)]=f(h,1-h)+(g(h,1)+g(0,1-h))/h2;
	for (int i=1;i<n-2;i++){
		b[(n-2)*(n-1)+i]=f(i*h,1-h)+g(i*h,1)/h2;
	}
	b[(n-2)*(n-1)+n-2]=f(1-h,1-h)+(g(1-h,1)+g(1,1-h))/h2;
	
	//Iterationsverfahren
	for (int k=0;k<iterNum;k++){
		for (int i=0;i<(n-1)*(n-1);i++){
			double dum=0; //das, was in der Klammer steht ("premature optimization"...)
			//aufgrund der Gestalt der Matrix vereinfachen sich die Summen drastisch
			dum=dum+4*u[i]/h2;
			if (i%(n-1)!=0){
				dum=dum-u[i-1]/h2;
			}
			if (i%(n-1)!=(n-2)){
				dum=dum-u[i+1]/h2;
			}
			if (i>n-2){
				dum=dum-u[i-n+1]/h2;
			}
			if (i<(n-1)*(n-2)){
				dum=dum-u[i+n-1]/h2;
			}
			dum=dum-b[i];
			u[i]=u[i]-omega*h2/4*dum;
		}
	}
	
	FILE* dat=fopen("./poisson.txt","w");
	//Ausgabe so, dass die Werte wie in einem mathematischen Koordinatensystem dastehen
	for (int i=n-2;i>-1;i--){
		for (int j=0;j<n-1;j++){
			fprintf(dat,"%f",u[i*(n-1)+j]);
			if (j<n-2){
				fprintf(dat,"\t");
			}
		}
		fprintf(dat,"\n");
	}
	fclose(dat);
	
	return 0;
}