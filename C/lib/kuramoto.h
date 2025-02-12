#ifndef KURAMATO_H
#define KURAMATO_H

// Kuramoto model

void kuramoto_euler	// Euler method
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_eulerpl // Euler method with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	const   double* const a,   // phase lags (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	const   double* const a,   // phase lags (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

#endif // KURAMATO_H
