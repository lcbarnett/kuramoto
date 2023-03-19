#ifndef MT64_H
#define MT64_H

// 64-bit Mersenne Twister PRNG (thread-safe)

#include <inttypes.h>

typedef uint64_t mtuint_t;

#define NN 312
typedef struct {
    mtuint_t mt[NN];
    int mti;
    int    iset;
    double gset;
} mt_t;
#undef NN

// initializes state
mtuint_t mt_seed(mt_t* const pstate, mtuint_t seed);

// copy full state of rng (for save/restore)
void mt_copy(mt_t* const dstate, const mt_t* const sstate);

// generates a random number on [0, 2^64-1]-interval
mtuint_t mt_uint(mt_t* const pstate);

// generates a random number uniform on [0,1)-real-interval
inline double mt_rand(mt_t* const pstate)
{
    return (double)(mt_uint(pstate) >> 11) * (1.0/9007199254740992.0);
}

// generates a random number uniform on (0,1)-real-interval
inline double mt_rand_pos(mt_t* const pstate)
{
	double x;
	do x = (double)(mt_uint(pstate) >> 11) * (1.0/9007199254740992.0); while (x == 0.0);
	return x;
}

// return integer in range 0,...,range-1 (be careful of overflow!!)
#define RANDI(itype,range,pstate) (itype)((double)(range)*mt_rand(pstate));

// generates a normally distributed random number from N(0,1)
double mt_randn(mt_t* const pstate);

// generates a Gamma variate from Gamma(a,b) [Marsaglia and Tsang]
double mt_rang(mt_t* const pstate, const double a, const double b);

// generates a standard Cauchy random variate
double mt_randc(mt_t* const pstate);

#endif // MT64_H
