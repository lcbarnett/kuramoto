#ifndef STULAN_H
#define STULAN_H

void stulan_euler // Euler method
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double        dt, // time integration step
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	const double* const a,  // dt*(growth constants) - (K summed over 2nd index)
	double*       const x,  // oscillator real part, to be computed by numerical ODE (pre-initialised with input)
	double*       const y   // oscillator imag part, to be computed by numerical ODE (pre-initialised with input)
);

#endif // STULAN_H
