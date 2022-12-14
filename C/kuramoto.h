#ifndef KURAMATO_H
#define KURAMATO_H

void kuramoto_euler	// Euler method (fast, less accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double        a,  // phase-lag (scalar)
	double*       const h   // oscillator phases computed by numerical ODE, pre-initialised
);

void kuramoto_rk4 // Classic Runge-Kutta ("RK4" - slower, more accurate)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*(coupling constants)
	const double        a,  // phase-lag (scalar)
	double*       const h   // oscillator phases computed by numerical ODE, pre-initialised
);

#endif // KURAMATO_H
