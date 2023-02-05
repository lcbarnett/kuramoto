#ifndef KUTILS_H
#define KUTILS_H

#include <stdint.h>
#include <error.h>
#include <stdlib.h>

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

// Uniform random (IEEE-754 53-bit resolution) double on [0,1)
//
// Note: uses POSIX random(), which is non-reantrant
// so not thread-safe/. Should probably use a better
// PRNG, like Mersenne Twister.

static inline double randu()
{
	// random() only returns 32 random bits - glue two together.
	return ((((uint64_t)random())|(((uint64_t)random())<<32))>>11)*(1.0/9007199254740992.0);
}

// Standard normal random double (Box-Muller, non-reantrant so not thread-safe)

double randn();

// get a random random seed

unsigned get_rand_seed();

// A basic stopwatch

double timer_start(const char mesg[]);
void   timer_stop(const double ts);

// Linear PCM encoding

int pcm_write(FILE* const fp, const double* const x, const size_t n, const int pcm, const double amin, const double amax);

#endif // KUTILS_H
