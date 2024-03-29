#ifndef KUTILS_H
#define KUTILS_H

#include <stdint.h>
#include <error.h>
#include <stdlib.h>
#include <limits.h>

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

// Linear PCM encoding

int pcm_write(FILE* const fp, const double* const x, const size_t n, const int pcm, const double amin, const double amax);

#endif // KUTILS_H
