#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h> // for malloc, etc.

// NOTE:  C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.

void kuramoto_euler	// Euler method
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
)
{
	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double hti = ht[i];
			double dhti = w[i];
			for (size_t j=0; j<N; ++j) dhti += Ki[j]*sin(ht[j]-hti);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_eulerpl // Euler method with phase lags
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double*       a,  // phase lags
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
)
{
	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double dhti = w[i];
			for (size_t j=0; j<N; ++j) dhti += Ki[j]*sin(ht[j]-hti-ai[j]);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const size_t        N, // number of oscillators
	const size_t        n, // number of integration increments
	const double* const w, // dt*frequencies
	const double* const K, // dt*(coupling constants)
	double*       const h, // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
	double*             k1 // buffer for RK4 coefficients (size must be 4*N)
)
{
	double* const k2 = k1+N;
	double* const k3 = k2+N;
	double* const k4 = k3+N;

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double hti = ht[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]-hti);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k1[j]-htk1i);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k2[j]-htk2i);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k3[j]-htk3i);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
		}
	}
}

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const size_t        N, // number of oscillators
	const size_t        n, // number of integration increments
	const double* const w, // dt*frequencies
	const double* const K, // dt*(coupling constants)
	const double*       a, // phase lags
	double*       const h, // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
	double*             k1 // buffer for RK4 coefficients (size must be 4*N)
)
{
	double* const k2 = k1+N;
	double* const k3 = k2+N;
	double* const k4 = k3+N;

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]-hti-ai[j]);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double* const ai = a+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k1[j]-htk1i-ai[j]);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double* const ai = a+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k2[j]-htk2i-ai[j]);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double* const ai = a+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(ht[j]+k3[j]-htk3i-ai[j]);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
		}
	}
}

void order_param // calculate order parameter magnitude
(
	const size_t N,
	const size_t n,
	const double* const h,
	double* const r
)
{
	const double OON = 1.0/(double)N;
	double* rt = r;
	for (const double* ht=h; ht<h+N*n; ht+=N,++rt) {
		double cmean = 0.0;
		double smean = 0.0;
		for (size_t i=0; i<N; ++i) {
#ifdef _GNU_SOURCE
			double c,s;
			sincos(ht[i],&s,&c);
			cmean += c;
			smean += s;
#else
			cmean += cos(ht[i]);
			smean += sin(ht[i]);
#endif
		}
		cmean *= OON;
		smean *= OON;
		*rt = hypot(cmean,smean);
	}
}

static inline double wmpi2pi(const double x) // wrap to [-pi,pi)
{
	return x > 0.0 ? fmod(x+M_PI,2.0*M_PI)-M_PI : fmod(x-M_PI,2.0*M_PI)+M_PI;
}

void phase_wrap // wrap all oscillator phases
(
	const size_t N,
	const size_t n,
	const double* const h,
	double* const theta
)
{
	for (size_t k=0; k<N*n; ++k) {
		theta[k] = wmpi2pi(h[k]);
	}
}
