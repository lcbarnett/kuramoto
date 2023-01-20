#ifndef KURAMATO_H
#define KURAMATO_H

void kuramoto_euler	// Euler method
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

void kuramoto_eulerpl // Euler method with phase lags
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double*       a,  // phase lags
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double*       a,  // phase lags
	double*       const h   // oscillator phases, to be computed by numerical ODE (pre-initialised with input)
);

#endif // KURAMATO_H
