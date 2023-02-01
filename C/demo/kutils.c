#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#ifdef __unix__
#include <unistd.h>
#include <time.h>
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

double timer_start(const char mesg[])
{
	printf("%s ...",mesg);
	fflush(stdout);
#ifdef __unix__
	return (double)clock()/(double)CLOCKS_PER_SEC;
#else
	return 0.0/0.0; // NaN
#endif
}

void timer_stop(const double ts)
{
#ifdef __unix__
	const double te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf(" %.4f seconds\n\n",te-ts);
#else
	printf(" done\n\n");
#endif
}

// Linear PCM (remember to free returned buffer after use!)

#define O16 ((uint16_t)1)
#define O32 ((uint32_t)1)
int pcm_write(FILE* const fp, const double* const x, const size_t n, const int pcm, const double amax, const double amin)
{
	if ((pcm != -32)&&(pcm != -64)&&(pcm != 16)&&(pcm != 24)) {
		fprintf(stderr,"ERROR: PCM must be unsigned 16- or 24-bit, or floating-point 32- or 64-bit\n");
		return 1;
	}
	if (pcm == -32) { // single-precision floating-point
		float* const u = calloc(n,sizeof(float)); // caller must free!!!
		float* uu = u;
		for (const double* xx=x; xx<x+n; ++xx) {
			*uu++ = (float)*xx;
		}
		if (fwrite(u,sizeof(float),n,fp) != n) {
			perror("ERROR: PCM write failed");
			free(u);
			return 2;
		}
		free(u);
	}
	else if  (pcm == -64) { // double-precision floating-point
		if (fwrite(x,sizeof(double),n,fp) != n) {
			perror("ERROR: PCM write failed");
			return 2;
		}
	}
	else {
		const size_t nbytes = (pcm/8)*n; // number of PCM bytes
		uchar_t* const u = calloc(nbytes,sizeof(uchar_t)); // caller must free!!!
		if (pcm == 16) {
			const uint16_t lomask = (O16<<8)-O16; // mask for low byte
			const double   maxval = (double)((O16<<16)-O16); // 16-bit max (double-precision floating-point)
			uchar_t* uu = u;
			for (const double* xx=x; xx<x+n; ++xx) {
				const uint16_t xpcm = (uint16_t)(maxval*((*xx-amin)/(amax-amin)));
				*uu++ = (uchar_t)((xpcm>>0)&lomask);
				*uu++ = (uchar_t)((xpcm>>8)&lomask);
			}
		}
		else { // pcm == 24
			const uint32_t lomask = (O32<<8)-O32; // mask for low byte
			const double   maxval = (double)((O32<<24)-O32); // 24-bit max (double-precision floating-point)
			uchar_t* uu = u;
			for (const double* xx=x; xx<x+n; ++xx) {
				const uint32_t xpcm = (uint32_t)(maxval*((*xx-amin)/(amax-amin)));
				*uu++ = (uchar_t)((xpcm>> 0)&lomask);
				*uu++ = (uchar_t)((xpcm>> 8)&lomask);
				*uu++ = (uchar_t)((xpcm>>16)&lomask);
			}
		}
		if (fwrite(u,sizeof(uchar_t),nbytes,fp) != nbytes) {
			perror("ERROR: PCM write failed");
			free(u);
			return 2;
		}
		free(u);
	}
	return 0;
}
#undef O32
#undef O16
