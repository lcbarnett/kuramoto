#include <stdio.h>
#include <math.h>

#include "kutils.h"

// Standard normal random double (Box-Muller, non-reantrant so not thread-safe)

double randn()
{
    static int    iset = 0;
    static double gset = 0.0;
    if (iset) {
	    iset=0;
	    return gset;
    }
    double v1,v2,rsq;
    do {
	    v1 = 2.0*randu()-1.0;
	    v2 = 2.0*randu()-1.0;
	    rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    const double fac = sqrt(-2.0*log(rsq)/rsq);
    gset = fac*v1;
    iset = 1;
    return fac*v2;
}

// get a random random seed (only implemented for Linux at the moment)

unsigned get_rand_seed()
{
	unsigned seed;
#ifdef	__linux__
	// if Linux, we seed from /dev/urandom
	if (getrandom(&seed,sizeof(unsigned),GRND_NONBLOCK) != sizeof(unsigned)) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
#else
	// else just set to 1
	seed = 1;
#endif
	return seed;
}
