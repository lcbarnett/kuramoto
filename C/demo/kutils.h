#ifndef KUTILS_H
#define KUTILS_H

// for random number generation from OS

#include <stdlib.h>
#ifdef	__linux__
#include <sys/random.h>
#endif

// Gnuplot default terminal

#ifdef _GNUPLOT_HAVE_PIPE
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


#endif // KUTILS_H
