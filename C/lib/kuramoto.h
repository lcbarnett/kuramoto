#ifndef KURAMATO_H
#define KURAMATO_H

#define TWOPI (2.0*M_PI)

// Kuramoto model

void kuramoto_euler	// Euler method
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (radians/sec)
	double* const  K,  // coupling constants (radians/sec)
	double* const  h   // oscillator phases (radians), initialised with input
);

void kuramoto_eulerpl // Euler method with phase lags
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (radians/sec)
	double* const  K,  // coupling constants (radians/sec)
	double* const  a,  // phase lags (radians)
	double* const  h   // oscillator phases (radians), initialised with input
);

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (radians/sec)
	double* const  K,  // coupling constants (radians/sec)
	double* const  h,  // oscillator phases (radians), initialised with input
	double* const  k1  // buffer for RK4 coefficients (size must be 4*N)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (radians/sec)
	double* const  K,  // coupling constants (radians/sec)
	double* const  a,  // phase lags (radians)
	double* const  h,  // oscillator phases (radians), initialised with input
	double* const  k1  // buffer for RK4 coefficients (size must be 4*N)
);

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const size_t N,        // number of oscillators
	const size_t n,        // number of integration increments
	const double* const h, // oscillator phases
	double* const r,       // order parameter magnitude
	double* const psi      // order parameter phase (NULL if not required)
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

// Utilities

static inline double wmpi2pi(const double x) // wrap to [-pi,pi)
{
	return x > 0.0 ? fmod(x+M_PI,2.0*M_PI)-M_PI : fmod(x-M_PI,2.0*M_PI)+M_PI;
}

void phase_wrap(const size_t m, double* const h);

#endif // KURAMATO_H
