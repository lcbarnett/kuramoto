#ifndef STULAN_H
#define STULAN_H

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

#endif // STULAN_H
