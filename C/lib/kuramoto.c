#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h>  // for malloc, etc.

#include "kuramoto.h"

// NOTE 1: C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.
//
// NOTE 2: All phases (including h) are dimensionless; multiply by 2π for radians.
//
// NOTE 3: To initialise input phase variables with Wiener noise, scale by sqrt(dt).
// For other (deterministic) input scale by dt.

void kuramoto_euler	// Euler method
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double hti = ht[i];
			double dhti = wdt[i];
			for (size_t j=0; j<N; ++j) dhti += Kdti[j]*sin(ht[j]-hti);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_eulerpl // Euler method with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	const   double* const a,   // phase lags (dimensionless)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double dhti = wdt[i];
			for (size_t j=0; j<N; ++j) dhti += Kdti[j]*sin(ht[j]-hti-ai[j]);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t  N,          // number of oscillators
	const   size_t  n,          // number of integration increments
	const   double* const  wdt, // frequencies x dt (dimensionless)
	const   double* const  Kdt, // coupling constants x dt (dimensionless)
	double* const   h           // oscillator phases, initialised with input (dimensionless)
)
{
	// coefficients buffers

	double k1[N],k2[N],k3[N],k4[N];

	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double hti = ht[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]-hti);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k1[j]-htk1i);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k2[j]-htk2i);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k3[j]-htk3i);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
	}
}

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	const   double* const a,   // phase lags (dimensionless)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	// coefficients buffers

	double k1[N],k2[N],k3[N],k4[N];

	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]-hti-ai[j]);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k1[j]-htk1i-ai[j]);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k2[j]-htk2i-ai[j]);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k3[j]-htk3i-ai[j]);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
	}
}
