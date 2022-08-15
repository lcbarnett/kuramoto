#include <math.h>   // for maths functions
#include <string.h> // for memcpy
#include <stdlib.h> // for malloc, etc.
#include <stdio.h>

void kuramoto_euler	// Euler method (fast, less accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double* const h0, // initial oscillator phases
	double*       const h   // oscillator phases computed by numerical ODE
)
{
	memcpy(h,h0,N*sizeof(double)); // initialise

	for (size_t t=0; t<n-1; ++t) {
		const double* const ht = h+N*t;
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			const double hti = ht[i];
			double ht1i = hti+w[i];
			for (size_t j=0; j<N; ++j) ht1i += K[j]*sin(ht[j]-hti);
			ht1[i] = ht1i; // update next time step
		}
	}
}

void kuramoto_rk4 // Classic Runge-Kutta ("RK4" - slower, more accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double* const h0, // initial oscillator phases
	double*       const h   // oscillator phases computed by numerical ODE
)
{
	// allocate buffer for intermediates (k1, k2, k3, k4)

	double* const kbuff = calloc(4*N,sizeof(double));
	double* const k1dt = kbuff;
	double* const k2dt = kbuff+N;
	double* const k3dt = kbuff+2*N;
	double* const k4dt = kbuff+3*N;

	// classic Runge-Kutta ("RK4" - slower, more accurate)

	memcpy(h,h0,N*sizeof(double)); // initialise

	for (size_t t=0; t<n-1; ++t) {
		const double* const ht = h+N*t;

		// k1
		for (size_t i=0; i<N; ++i) {
			const double hti = ht[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += K[j]*sin(ht[j]-hti);
			k1dt[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double htpk1dti = ht[i]+k1dt[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += K[j]*sin(ht[j]+k1dt[j]-htpk1dti);
			k2dt[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double htpk2dti = ht[i]+k2dt[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += K[j]*sin(ht[j]+k2dt[j]-htpk2dti);
			k3dt[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double htpk3dti = ht[i]+k3dt[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += K[j]*sin(ht[j]+k3dt[j]-htpk3dti);
			k4dt[i] = ki;
		}

		// update next time step
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			ht1[i] = ht[i] + (k1dt[i]+4.0*k2dt[i]+4.0*k3dt[i]+k4dt[i])/6.0;
		}
	}

	// free buffer

	free(kbuff);
}
