#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.

// NOTE:  C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrix K before calling.

void kuramoto_euler	// Euler method (fast, less accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double        a,  // phase-lag (scalar)
	double*       const h   // oscillator phases computed by numerical ODE, pre-initialised
)
{
	for (size_t t=0; t<n-1; ++t) {
		const double* const ht = h+N*t;
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htipa = ht[i]+a;
			double ht1i = ht[i]+w[i];
			for (size_t j=0; j<N; ++j) ht1i += Ki[j]*sin(ht[j]-htipa);
			ht1[i] += ht1i; // update next time step (adding in input already in ht1)
		}
	}

}

void kuramoto_rk4 // Classic Runge-Kutta ("RK4" - slower, more accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double        a,  // phase-lag (scalar)
	double*       const h   // oscillator phases computed by numerical ODE, pre-initialised
)
{
	// allocate buffer for intermediates (k1, k2, k3, k4)

	double* const kbuff = calloc(4*N,sizeof(double));
	double* const k1dt = kbuff;
	double* const k2dt = kbuff+N;
	double* const k3dt = kbuff+2*N;
	double* const k4dt = kbuff+3*N;

	for (size_t t=0; t<n-1; ++t) {
		const double* const ht = h+N*t;

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htipa = ht[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]-htipa);
			k1dt[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htpk1dtipa = ht[i]+k1dt[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k1dt[j]-htpk1dtipa);
			k2dt[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htpk2dtipa = ht[i]+k2dt[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k2dt[j]-htpk2dtipa);
			k3dt[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htpk3dtipa = ht[i]+k3dt[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k3dt[j]-htpk3dtipa);
			k4dt[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			ht1[i] += ht[i] + (k1dt[i]+4.0*k2dt[i]+4.0*k3dt[i]+k4dt[i])/6.0;
		}
	}

	// free buffer

	free(kbuff);
}
