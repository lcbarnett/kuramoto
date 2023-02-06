#ifndef KURAMATO_H
#define KURAMATO_H

#define TWOPI (2.0*M_PI)

void kuramoto_euler	// Euler method
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

void kuramoto_eulerpl // Euler method with phase lags
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	const double* const a,  // phase lags
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const size_t        N, // number of oscillators
	const size_t        n, // number of integration increments
	const double* const w, // dt*frequencies
	const double* const K, // dt*frequencies*(coupling constants)/N
	double*       const h, // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
	double*       const k1 // buffer for RK4 coefficients (size must be 4*N)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const size_t        N, // number of oscillators
	const size_t        n, // number of integration increments
	const double* const w, // dt*frequencies
	const double* const K, // dt*frequencies*(coupling constants)/N
	const double* const a, // phase lags
	double*       const h, // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
	double*       const k1 // buffer for RK4 coefficients (size must be 4*N)
);

void order_param // calculate order parameter magnitude
(
	const size_t N,        // number of oscillators
	const size_t n,        // number of integration increments
	const double* const h, // oscillator phases
	double* const r,       // order parameter magnitude
	double* const psi      // order parameter phase (NULL if not required)
);

void phase_wrap(const size_t m, double* const h);

#endif // KURAMATO_H
