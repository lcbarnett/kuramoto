#ifndef KURAMATO_H
#define KURAMATO_H

#define TWOPI (2.0*M_PI)

// Stuart-Landau: expected amplitude of sum of N oscillators
// uniform random on unit circle is approx SLMAGIC*sqrt(N)

#define SLMAGIC (0.887)

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
	double* const         k1,  // buffer for RK4 coefficients (size must be 4*N)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	const   double* const a,   // phase lags (radians)
	double* const         k1,  // buffer for RK4 coefficients (size must be 4*N)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
);

// Stuart-Landau model

void stulan_euler // Euler method
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (1/sec)
	double* const  K,  // coupling constants (1/sec)
	double* const  a,  // growth constants (1/sec)
	double* const  x,  // oscillator real part (dimensionless), initialised with input
	double* const  y   // oscillator imag part (dimensionless), initialised with input
);

void stulan_magnitudes // calculate magnitudes of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r  // oscillator magnitude
);

void stulan_phases // calculate phases of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  h  // oscillator phase
);

void stulan_order_param // calculate order parameter magnitude/phase
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r, // order parameter magnitude
	double* const  p  // order parameter phase (NULL if not required)
);

// Utilities

static inline double phasewrap(const double x, const double u) // wrap to [-u,u)
{
	return x > 0.0 ? fmod(x+u,2.0*u)-u : fmod(x-u,2.0*u)+u;
}

void phase_wrap(const size_t m, double* const h, const double u); // wrap to [-u,u)

#endif // KURAMATO_H
