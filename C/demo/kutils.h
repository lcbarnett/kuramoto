#ifndef KUTILS_H
#define KUTILS_H

#include <stdint.h>

// for random number generation from OS

#include <stdlib.h>
#ifdef	__linux__
#include <sys/random.h>
#endif

// Gnuplot default terminal

#ifdef _HAVE_GNUPLOT
#ifdef __unix__
#define GPTERM "x11"
#endif
#ifdef __APPLE__
#define GPTERM "aqua"
#endif
#endif

// Uniform random double on [0,1) [Note: you might want a better PRNG :-)]

static inline double randu()
{
	return (double)random()/((double)(RAND_MAX)+1.0); // (random() is non-reantrant so not thread-safe)
}

// Standard normal random double (Box-Muller, non-reantrant so not thread-safe)

double randn();

// get a random random seed (only implemented for Linux at the moment)

unsigned get_rand_seed();

// A basic stopwatch

double timer_start(const char mesg[]);
void   timer_stop(const double ts);

// Linear PCM

static inline uint16_t pcm16(const double x, const double amax, const double amin)
{
	return (uint16_t)(((double)((((uint16_t)1)<<16)-1))*((x-amin)/(amax-amin)));
}

static inline uint32_t pcm24(const double x, const double amax, const double amin)
{
	return (uint32_t)(((double)((((uint32_t)1)<<24)-1))*((x-amin)/(amax-amin)));
}

void xpcm16(const double* const x, uint16_t* const u, const size_t n, const double amax, const double amin);
void xpcm24(const double* const x, uint32_t* const u, const size_t n, const double amax, const double amin);

#endif // KUTILS_H
