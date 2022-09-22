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
	for (const double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double hti = ht[i];
			const double htai = hti+a;
			double ht1i = hti+w[i];
			for (size_t j=0; j<N; ++j) ht1i += Ki[j]*sin(ht[j]-htai);
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
	double* const k1 = kbuff;
	double* const k2 = kbuff+N;
	double* const k3 = kbuff+2*N;
	double* const k4 = kbuff+3*N;

	for (const double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htai = ht[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]-htai);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk1ai = ht[i]+k1[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k1[j]-htk1ai);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk2ai = ht[i]+k2[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k2[j]-htk2ai);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk3ai = ht[i]+k3[i]+a;
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k3[j]-htk3ai);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = (double* const)ht+N;
		for (size_t i=0; i<N; ++i) {
			ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
		}
	}

	// free buffer

	free(kbuff);
}
