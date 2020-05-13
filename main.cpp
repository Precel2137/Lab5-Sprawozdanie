#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "winbgi2.h"
#include "rk4.h"
#define l 0.4
#define g 9.81


void rhs_fun(double t, double *X, double *F)
{
  F[0] = X[1];
  F[1] = -sin(X[0])*g/l;


}
void veuler(double t, double *X, double h, int n, void (* fun)(double, double *,double *), double *X1)
{
	double F[2];
	fun(t,X,F);
	X1[0] = X[0] + h*F[0];	
	X1[1] = X[1] + h*F[1];

}

void main()
{
	graphics(800, 800);
	scale(-180, -10, 180, 10);
	double alfa,omega,x[2],x1[2],t0=0.0,tk=10.0,h=0.005,m = 2.5,E,Ek,Ep;
	printf("Podaj wartosc poczatkowa alfa(w radianach):\n");
	scanf("%lf",&x[0]);
	printf("Podaj wartosc poczatkowa omega\n");
	scanf("%lf",&x[1]);
	alfa = x[0];
	omega = x[1];
	x1[0] = x[0];
	x1[1] = x[1];
	while(t0 <= tk)
	{
		setcolor(1.);
		point(x1[0]*180./M_PI,x1[1]);
		//veuler(t0,x,h,2,rhs_fun,x1);
		vrk4(t0,x,h,2,rhs_fun,x1);
		x[0] = x1[0];
		x[1] = x1[1];
		t0 += h;


	}
	t0 = 0.;
	printf("Wpisz dowolna wartosc i wcisnij enter aby wyswietlic wykres rozkladu energii w czasie\n");
	scanf("%lf",&x[0]);
	x[0]=alfa;
	x[1]=omega;
	setgray(1.);
	scale(0, 0, 10, 10);
	while(t0 <= tk)
	{
		setcolor(0.);
		point(t0,0);
		//veuler(t0,x,h,2,rhs_fun,x1);
		vrk4(t0,x,h,2,rhs_fun,x1);
		Ek = m*pow(l,2.)/2*pow(x[1],2.);
		Ep = m*g*l*(1-cos(x[0]));
		E = Ep + Ek;
		setcolor(0.5);
		point(t0,Ek);
		setcolor(1.);
		point(t0,Ep);
		setcolor(0.75);
		point(t0,E);
		x[0] = x1[0];
		x[1] = x1[1];
		t0 += h;


	}
	
	
	
	
	
	printf("program zakonczyl dzialanie");
	scanf("%lf",&x[0]);
	wait();
}

