#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h>  // for printf, etc.

#include "kuramoto.h"

// NOTE 1: C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.
//
// NOTE 2: All phases (including h) are dimensionless; multiply by 2π for radians.
//
// NOTE 3: To initialise input phase variables with Wiener noise, scale by sqrt(dt).
// For other (deterministic) input scale by dt.

void kmoto_euler
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	double* const         x  // oscillator phases, initialised with input
)
{
	printf("\nkmoto_euler\n");
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		kmoto_fun(udot,u,N,w,K);
		double* const unext = u+N;
		for (size_t i=0; i<N; ++i) unext[i] += u[i] + h*udot[i];
	}
}

void kmoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	double* const         x  // oscillator phases, initialised with input
)
{
	printf("\nkmoto_rk4\n");
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		kmoto_fun(udot1,u,N,w,K);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		kmoto_fun(udot2,v,N,w,K);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		kmoto_fun(udot3,v,N,w,K);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		kmoto_fun(udot4,v,N,w,K);
		double* const unext = u+N;
		for (size_t i=0; i<N; ++i) unext[i] += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

void kmotopl_euler
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	const   double* const a, // phase lags
	double* const         x  // oscillator phases, initialised with input
)
{
	printf("\nkmotopl_euler\n");
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		kmotopl_fun(udot,u,N,w,K,a);
		double* const unext = u+N;
		for (size_t i=0; i<N; ++i) unext[i] += u[i] + h*udot[i];
	}
}

void kmotopl_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	const   double* const a, // phase lags
	double* const         x  // oscillator phases, initialised with input
)
{
	printf("\nkmotopl_rk4\n");
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		kmotopl_fun(udot1,u,N,w,K,a);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		kmotopl_fun(udot2,v,N,w,K,a);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		kmotopl_fun(udot3,v,N,w,K,a);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		kmotopl_fun(udot4,v,N,w,K,a);
		double* const unext = u+N;
		for (size_t i=0; i<N; ++i) unext[i] += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

void rossler_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nrossler_euler\n");
	const size_t N = 3;
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		rossler_fun(udot,u,a,b,c);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i];
	}
}

void rossler_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nrossler_heun\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	double udot1[N], udot2[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		rossler_fun(udot1,u,a,b,c);
		for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i];
		rossler_fun(udot2,v,a,b,c);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]);
	}
}

void rossler_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nrossler_rk4\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		rossler_fun(udot1,u,a,b,c);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		rossler_fun(udot2,v,a,b,c);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		rossler_fun(udot3,v,a,b,c);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		rossler_fun(udot4,v,a,b,c);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

void lorenz_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlorenz_euler\n");
	const size_t N = 3;
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lorenz_fun(udot,u,s,r,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i];
	}
}

void lorenz_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlorenz_heun\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	double udot1[N], udot2[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lorenz_fun(udot1,u,s,r,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i];
		lorenz_fun(udot2,v,s,r,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]);
	}
}

void lorenz_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlorenz_rk4\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lorenz_fun(udot1,u,s,r,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		lorenz_fun(udot2,v,s,r,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		lorenz_fun(udot3,v,s,r,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		lorenz_fun(udot4,v,s,r,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

void thomas_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nthomas_euler\n");
	const size_t N = 3;
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		thomas_fun(udot,u,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i];
	}
}

void thomas_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nthomas_heun\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	double udot1[N], udot2[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		thomas_fun(udot1,u,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i];
		thomas_fun(udot2,v,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]);
	}
}

void thomas_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nthomas_rk4\n");
	const size_t N = 3;
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		thomas_fun(udot1,u,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		thomas_fun(udot2,v,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		thomas_fun(udot3,v,b);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		thomas_fun(udot4,v,b);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

void lrnz96_euler
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlrnz96_euler\n");
	double udot[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lrnz96_fun(udot,u,N,F);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i];
	}
}

void lrnz96_heun
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlrnz96_heun\n");
	const double h2 = h/2.0;
	double udot1[N], udot2[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lrnz96_fun(udot1,u,N,F);
		for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i];
		lrnz96_fun(udot2,v,N,F);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]);
	}
}

void lrnz96_rk4
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
)
{
	printf("\nlrnz96_rk4\n");
	const double h2 = h/2.0;
	const double h6 = h/6.0;
	double udot1[N],udot2[N],udot3[N],udot4[N];
	double v[N];
	for (double* u=x; u<x+N*(n-1); u+=N) {
		lrnz96_fun(udot1,u,N,F);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i];
		lrnz96_fun(udot2,v,N,F);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i];
		lrnz96_fun(udot3,v,N,F);
		for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i];
		lrnz96_fun(udot4,v,N,F);
		double* const u1 = u+N;
		for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]);
	}
}

// Utilities

void phase_wrap(const size_t m, double* const h, const double u) // wrap to [-u,u)
{
	for (size_t k=0; k<m; ++k) h[k] = phasewrap(h[k],u);
}

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
)
{
	const double OON = 1.0/(double)N;
	double* rt = r;
	double* pt = p;
	for (const double* ht=h; ht<h+N*n; ht+=N) {
		double xmt = 0.0;
		double ymt = 0.0;
		for (size_t i=0; i<N; ++i) {
#ifdef _GNU_SOURCE
			double c,s;
			sincos(ht[i],&s,&c);
			xmt += c;
			ymt += s;
#else
			xmt += cos(ht[i]);
			ymt += sin(ht[i]);
#endif
		}
		xmt *= OON;
		ymt *= OON;
		*rt++ = hypot(xmt,ymt);
		if (p) *pt++ = atan2(ymt,xmt);
	}
}
