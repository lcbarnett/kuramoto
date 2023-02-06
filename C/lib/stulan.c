#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h> // for malloc, etc.

// NOTE:  C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.

void stulan_euler // Euler method
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double        dt, // time integration step
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	const double* const a,  // dt*(growth constants) - (K summed over 2nd index)
	double*       const x,  // oscillator real part, to be computed by numerical ODE (pre-initialised with input)
	double*       const y   // oscillator imag part, to be computed by numerical ODE (pre-initialised with input)
)
{
	double* yt=y;
	for (double* xt=x; xt<x+N*(n-1); xt+=N,yt+=N) {
		double* const xt1 = xt+N;
		double* const yt1 = yt+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xt[i]*xt[i]+yt[i]*yt[i]);
			double dxti = vi*xti - w[i]*yt[i];
			double dyti = vi*yti + w[i]*xt[i];
			for (size_t j=0; j<N; ++j) dxti += Ki[j]*xt[j];
			for (size_t j=0; j<N; ++j) dyti += Ki[j]*yt[j];
			yt1[i] += yti+dyti; // update next time step (adding in input already in xt1)
			xt1[i] += xti+dxti; // update next time step (adding in input already in yt1)
		}
	}
}

/*
void stulan_rk4 // Classic Runge-Kutta (RK4)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double        dt, // time integration step
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	const double* const a,  // dt*(growth constants) - (K summed over 2nd index)
	double*       const x,  // oscillator real part, to be computed by numerical ODE (pre-initialised with input)
	double*       const y,  // oscillator imag part, to be computed by numerical ODE (pre-initialised with input)
	double*       const k   // buffer for RK4 coefficients (size must be 8*N)
)
{
	double* const kx1 = k  +N;
	double* const ky1 = kx1+N;
	double* const kx2 = ky1+N;
	double* const ky2 = kx2+N;
	double* const kx3 = ky2+N;
	double* const ky3 = kx3+N;
	double* const kx4 = ky3+N;
	double* const ky4 = kx4+N;

	double* yt=y;
	for (double* xt=x; xt<x+N*(n-1); xt+=N,yt+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xt[i]*xt[i]+yt[i]*yt[i]);
			double kxi = vi*xti - w[i]*yt[i];
			double kyi = vi*yti + w[i]*xt[i];
			for (size_t j=0; j<N; ++j) kxi += Ki[j]*xt[j];
			for (size_t j=0; j<N; ++j) kyi += Ki[j]*yt[j];
			kx1[i] = kxi;
			ky1[i] = kyi;

		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xt[i]*xt[i]+yt[i]*yt[i]);

			const double xtk1i = xt[i]+kx1[i];
			double kxi = vi*xti - w[i]*yt[i];
			double kyi = vi*yti + w[i]*xt[i];
			for (size_t j=0; j<N; ++j) kxi += Ki[j]*(xt[j]+kx1[j]-xtk1i);
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k1[j]-htk1i);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk2i = xt[i]+k2[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k2[j]-htk2i);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk3i = xt[i]+k3[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k3[j]-htk3i);
			k4[i] = ki;
		}

		// update next time step (adding in input already in xt1)
		double* const xt1 = xt+N;
		for (size_t i=0; i<N; ++i) {
			xt1[i] += xt[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
		}
	}
}
*/
