#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mt64.h"

#define NN 312
#define MM 156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM UINT64_C(0xFFFFFFFF80000000) // Most  significant 33 bits
#define LM UINT64_C(0x7FFFFFFF)         // Least significant 31 bits

// initializes state
mtuint_t mt_seed(mt_t* const pstate, mtuint_t seed)
{
//	You MUST seed the rng before it is used!!!

	if (seed == 0) { // initialise from /dev/urandom
		FILE* fp = fopen("/dev/urandom","r");
		if (fp == NULL)                              {perror("mt_seed: failed to open /dev/urandom" ); exit(EXIT_FAILURE);}
		if (fread(&seed,sizeof(mtuint_t),1,fp) != 1) {perror("mt_seed: failed to read /dev/urandom" ); exit(EXIT_FAILURE);}
		if (fclose(fp) != 0)                         {perror("mt_seed: failed to close /dev/urandom"); exit(EXIT_FAILURE);}
	}

    pstate->mt[0] = seed;
    for (pstate->mti=1; pstate->mti<NN; pstate->mti++) pstate->mt[pstate->mti] = (UINT64_C(6364136223846793005) * (pstate->mt[pstate->mti-1] ^ (pstate->mt[pstate->mti-1] >> 62)) + (mtuint_t)pstate->mti);

	pstate->iset = 0;
	pstate->gset = 0.0;

    return seed;
}

// deep copy of rng state (for save/restore)
void mt_copy(mt_t* const dstate, const mt_t* const sstate)
{
	memcpy(dstate->mt,sstate->mt,NN*sizeof(mtuint_t));
	dstate->mti  = sstate->mti;
	dstate->iset = sstate->iset;
	dstate->gset = sstate->gset;
}

// generates a random number on [0, 2^64-1]-interval
mtuint_t mt_uint(mt_t* const pstate)
{
    int i;
    mtuint_t x;
    static const mtuint_t mag01[2]={UINT64_C(0), MATRIX_A};

    if (pstate->mti >= NN) { // generate NN words at one time

        // if init_genrand64() has not been called, a default initial seed is used
        // if (pstate->mti == NN+1) mt64_seed(pstate,UINT64_C(5489));
        // No! We don't want this. The user MUST seed first! (probably segfault if they forget)

        for (i=0;i<NN-MM;i++) {
            x = (pstate->mt[i]&UM)|(pstate->mt[i+1]&LM);
            pstate->mt[i] = pstate->mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        for (;i<NN-1;i++) {
            x = (pstate->mt[i]&UM)|(pstate->mt[i+1]&LM);
            pstate->mt[i] = pstate->mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        x = (pstate->mt[NN-1]&UM)|(pstate->mt[0]&LM);
        pstate->mt[NN-1] = pstate->mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];

        pstate->mti = 0;
    }

    x = pstate->mt[pstate->mti++];

    x ^= (x >> 29) & UINT64_C(0x5555555555555555);
    x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
    x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
    x ^= (x >> 43);

    return x;
}

// generates a Gaussian variate from N(0,1)
double mt_randn(mt_t* const pstate)
{
    if (pstate->iset) {
	    pstate->iset=0;
	    return pstate->gset;
    }

    double v1,v2,rsq;
    do {
	    v1 = 2.0*mt_rand(pstate)-1.0;
	    v2 = 2.0*mt_rand(pstate)-1.0;
	    rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    const double fac = sqrt(-2.0*log(rsq)/rsq);
    pstate->gset = fac*v1;
    pstate->iset = 1;
    return fac*v2;
}

// generates a Gamma variate from Gamma(a,b) [Marsaglia and Tsang]
double mt_rang(mt_t* const pstate, const double a, const double b)
{
	if (a < 1) {
		double u = mt_rand_pos(pstate);
		return mt_rang(pstate,1.0+a,b)*pow (u,1.0/a);
	}

	const double d = a - 1.0/3.0;
	const double c = (1.0/3.0)/sqrt(d);
	double x,v,u;

	while (1) {
		do {
			x = mt_randn(pstate);
			v = 1.0 + c*x;
		} while (v <= 0);

		v = v*v*v;
		u = mt_rand_pos(pstate);

		if (u < 1.0-0.0331*x*x*x*x) break;

		if (log (u) < 0.5*x*x+d*(1.0-v+log (v))) break;
	}

	return (d*v)/b; // NOTE: scaling like Wikipedia - in gsl this would be b*d*v
}

// generates a standard Cauchy random variate
double mt_randc(mt_t* const pstate)
{
	 return tan(M_PI*(mt_rand_pos(pstate)-0.5));
}

#undef NN
#undef MM
#undef MATRIX_A
#undef UM
#undef LM
