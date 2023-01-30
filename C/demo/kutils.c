#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#ifdef __unix__
#include <unistd.h>
#endif

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

// get a random random seed (only implemented for Linux/POSIX at the moment)

unsigned get_rand_seed()
{
	unsigned seed;
#ifdef __unix__
	#ifdef	__linux__
	if (getrandom(&seed,sizeof(unsigned),GRND_NONBLOCK) != sizeof(unsigned)) {
		perror("random seed generation failed");
		return EXIT_FAILURE;
	}
	#else
	const int fd = open("/dev/urandom",O_RDONLY);
	if (fd < 0) {
		perror("random seed generation: failed to open /dev/urandom");
		return EXIT_FAILURE;
	}
	if (read(fd,&seed,sizeof(unsigned)) != sizeof(unsigned)) {
		perror("random seed generation: failed to read /dev/urandom");
		return EXIT_FAILURE;
	}
	if (close(fd) != 0) {
		perror("random seed generation: failed to close /dev/urandom");
		return EXIT_FAILURE;
	}
	#endif // __linux__
#else
	fprintf(stderr,"WARNING: no OS random seed - setting to 1\n");
	seed = 1;
#endif // __unix__
	return seed;
}
