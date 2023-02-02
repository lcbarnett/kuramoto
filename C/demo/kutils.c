#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#ifdef	__linux__
#include <sys/random.h> // for getrandom()
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
#endif
	return seed;
}

double timer_start(const char mesg[])
{
	printf("%s ...",mesg);
	fflush(stdout);
	return (double)clock()/(double)CLOCKS_PER_SEC;
}

void timer_stop(const double ts)
{
	const double te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf(" %.4f seconds\n\n",te-ts);
}

// Linear PCM (remember to free returned buffer after use!)

#define O16 ((uint16_t)1)
#define O32 ((uint32_t)1)
int pcm_write(FILE* const fp, const double* const x, const size_t n, const int pcm, const double amax, const double amin)
{
	if (pcm == 16) { // 2 bytes = uint16_t
		const double maxfac = (double)((O16<<16)-O16)/(amax-amin); // 16-bit max factor
		uint16_t* const u = calloc(n,sizeof(uint16_t));
		if (u == NULL) {
			perror("PCM memory allocation failed");
			free(u);
			return -1;
		}
		uint16_t* uu = u;
		for (const double* xx=x; xx<x+n; ++xx) {
			*uu++ = (uint16_t)(maxfac*(*xx-amin));
		}
		if (fwrite(u,sizeof(uint16_t),n,fp) != n) {
			perror("PCM write failed");
			free(u);
			return -1;
		}
		free(u);
		return 0;
	}

	if (pcm == 24) { // 3 bytes, but no native uint24_t type, so we use 3 unsigned chars
		const double maxfac = (double)((O32<<24)-O32)/(amax-amin); // 24-bit max factor
		const size_t nbytes = 3*n; // number of PCM bytes
		uchar_t* const u = calloc(nbytes,sizeof(uchar_t));
		if (u == NULL) {
			perror("PCM memory allocation failed");
			free(u);
			return -1;
		}
		const uint32_t lomask = (O32<<8)-O32; // mask for low byte
		uchar_t* uu = u;
		for (const double* xx=x; xx<x+n; ++xx) {
			const uint32_t xpcm = (uint32_t)(maxfac*(*xx-amin));
			*uu++ = (uchar_t)((xpcm>> 0)&lomask); // don't really need the shift
			*uu++ = (uchar_t)((xpcm>> 8)&lomask);
			*uu++ = (uchar_t)((xpcm>>16)&lomask); // don't really need the mask
		}
		if (fwrite(u,sizeof(uchar_t),nbytes,fp) != nbytes) {
			perror("PCM write failed");
			free(u);
			return -1;
		}
		free(u);
		return 0;
	}

	if (pcm == -32) { // single-precision floating-point
		float* const u = calloc(n,sizeof(float));
		if (u == NULL) {
			perror("PCM memory allocation failed");
			free(u);
			return -1;
		}
		float* uu = u;
		for (const double* xx=x; xx<x+n; ++xx) {
			*uu++ = (float)*xx;
		}
		if (fwrite(u,sizeof(float),n,fp) != n) {
			perror("PCM write failed");
			free(u);
			return -1;
		}
		free(u);
		return 0;
	}

	if  (pcm == -64) { // double-precision floating-point - got that already!
		if (fwrite(x,sizeof(double),n,fp) != n) {
			perror("PCM write failed");
			return -1;
		}
		return 0;
	}

	error(0,0,"PCM must be unsigned 16- or 24-bit, or floating-point 32- or 64-bit");
	return -1;
}
#undef O32
#undef O16
