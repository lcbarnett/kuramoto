#ifndef KUTILS_H
#define KUTILS_H

#include <stdint.h>
#include <error.h>
#include <stdlib.h>
#include <limits.h>

#define TWOPI (2.0*M_PI)

// Stuart-Landau: expected amplitude of sum of N oscillators
// uniform random on unit circle is approx SLMAGIC*sqrt(N)

#define SLMAGIC (0.887)

// Check that uint64_t is implemented

#if !defined(UINT64_MAX)
#error "Sorry, this code requires 64-bit unsigned integers"
#endif

// Gnuplot default terminal

#ifdef _HAVE_GNUPLOT
	#if defined(__unix__)
		#define GPTERM "x11"
	#elif defined(__APPLE__)
		#define GPTERM "aqua"
	#elif defined(_WIN32) ||  defined(_WIN64)
		#define GPTERM "windows"
	#else
		#error "Don't know how to set Gnuplot terminal"
	#endif
#endif

typedef unsigned char uchar_t;

/* random() and friends deprecated in favour of (thread-safe) Mersenne Twister

// Uniform random (IEEE-754 53-bit resolution) double on [0,1)
//
// Note: uses POSIX random(), which is non-reantrant
// so not thread-safe/. Should probably use a better
// PRNG, like Mersenne Twister.

static inline double randu()
{
	return (double)random()/(double)RAND_MAX;
}

// Uniform random (IEEE-754 53-bit resolution) double on (0,1)

static inline double randu0()
{
	double x;
	do x = randu(); while (x == 0.0);
	return x;
}

// Standard normal random double (Box-Muller, non-reantrant so not thread-safe)

double randn();

// Standard Cauchy random double

double randc();

// get a random random seed

unsigned get_rand_seed();

*/

// A basic stopwatch

double timer_start(const char mesg[]);
void   timer_stop(const double ts);

// Phase stuff

static inline double phasewrap(const double x, const double u) // wrap to [-u,u)
{
	return x > 0.0 ? fmod(x+u,2.0*u)-u : fmod(x-u,2.0*u)+u;
}

void phase_wrap(const size_t m, double* const h, const double u); // wrap to [-u,u)

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
);

// Linear PCM encoding

int pcm_write(FILE* const fp, const double* const x, const size_t n, const int pcm, const double amin, const double amax);

#endif // KUTILS_H
