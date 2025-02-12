#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h>  // for malloc, etc.

#include "stulan.h"

// NOTE 1: C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.
//
// NOTE 2: All phases (including h) are dimensionless; multiply by 2π for radians.
//
// NOTE 3: To initialise input phase variables with Wiener noise, scale by sqrt(dt).
// For other (deterministic) input scale by dt.


// Stuart-Landau model

void stulan_euler // Euler method
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (1/sec)
	double* const  K,  // coupling constants (1/sec)
	double* const  a,  // growth constants (1/sec)
	double* const  x,  // oscillator real part (dimensionless), initialised with input
	double* const  y   // oscillator imag part (dimensionless), initialised with input
)
{
	// scale parameters appropriately according to time increment

	const double sqrtdt = sqrt(dt);
	for (size_t i=0; i<N;   ++i) w[i] *= dt;
	for (size_t i=0; i<N*N; ++i) K[i] *= dt;
	for (size_t i=0; i<N;   ++i) a[i] *= dt;
	for (size_t i=N; i<N*n; ++i) x[i] *= sqrtdt; // cf. Ornstein-Uhlenbeck process
	for (size_t i=N; i<N*n; ++i) y[i] *= sqrtdt; // cf. Ornstein-Uhlenbeck process

	// adjust growth constants by mean coupling

	for (size_t i=0; i<N; ++i) {
		const double* const Ki = K+N*i;
		double Kbari = 0.0;
		for (size_t j=0; j<N; ++j) Kbari += Ki[j];
		a[i] -= Kbari;
	}

	// ODE solver

	for (double *xt=x, *yt = y; xt<x+N*(n-1); xt+=N, yt+=N) {
		double* const xt1 = xt+N;
		double* const yt1 = yt+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xti*xti+yti*yti);
			double dxti = vi*xti - w[i]*yti;
			double dyti = vi*yti + w[i]*xti;
			for (size_t j=0; j<N; ++j) dxti += Ki[j]*xt[j];
			for (size_t j=0; j<N; ++j) dyti += Ki[j]*yt[j];
			xt1[i] += xti+dxti; // update next time step (adding in input already in yt1)
			yt1[i] += yti+dyti; // update next time step (adding in input already in xt1)
		}
	}
}

void stulan_magnitudes // calculate magnitudes of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r  // oscillator magnitude
)
{
	double* rt = r;
	for (double *xt=x, *yt = y; xt<x+N*n; ++xt, ++yt) *rt++ = hypot(*xt,*yt);
}

void stulan_phases // calculate phases of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  h  // oscillator phase
)
{
	double* ht = h;
	for (double *xt=x, *yt = y; xt<x+N*n; ++xt, ++yt) *ht++ = atan2(*yt,*xt);
}

void stulan_order_param // calculate order parameter magnitude/phase
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r, // order parameter magnitude
	double* const  p  // order parameter phase (NULL if not required)
)
{
	const double OON = 1.0/(double)N;
	double* rt = r;
	double* pt = p;
	for (double *xt=x, *yt = y; xt<x+N*n; xt+=N, yt+=N) {
		double xmt = 0.0;
		for (size_t i=0; i<N; ++i)  xmt += xt[i];
		xmt *= OON;
		double ymt = 0.0;
		for (size_t i=0; i<N; ++i)  ymt += yt[i];
		ymt *= OON;
		*rt++ = hypot(xmt,ymt);
		if (p) *pt++ = atan2(ymt,xmt);
	}
}
