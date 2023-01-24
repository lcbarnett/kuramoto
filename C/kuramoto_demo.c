#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef	__linux__
#include <sys/random.h>
#endif
#include "kuramoto.h"

// Uniform random double on [0,1) [Note: you might want a better PRNG :-)]

static inline double randu()
{
	return (double)random()/((double)(RAND_MAX)+1.0); // (random() is non-reantrant so not thread-safe)
}

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

// Main function. Call as:
//
// kuramoto_demo <seed>

int main(int argc, char *argv[])
{
	// seed random number generator (from command line if you want predictability)

	const unsigned seed = argc > 1 ? (unsigned)atoi(argv[1]) : get_rand_seed();
	printf("\nrandom seed = %u\n\n",seed);
	srand(seed);

	// Kuramoto model size parameters

	const size_t N  = 4;                          // number of oscillators
	const double T  = 200.0;                      // total integration time
	const double dt = 0.01;                       // integration step size
	const size_t n  = (size_t)ceil(T/dt);         // number of integration steps

	// allocate memory

	double* const w = calloc(N,  sizeof(double)); // oscillator frequencies
	double* const K = calloc(N*N,sizeof(double)); // coupling constants
	double* const h = calloc(N*n,sizeof(double)); // oscillator phases, unwrapped
	double* const r = calloc(n,  sizeof(double)); // order parameter

	// set up some random frequencies

	const double wmean = 0.0;
	const double wsdev = M_PI/7.0;
	for (size_t i=0; i<N; ++i) {
		w[i] = dt*(wmean+wsdev*randn()); // scale frequencies by dt
	}

	// set up some random coupling constants

	const double Kmean = 0.8/(double)N;
	const double Ksdev = Kmean/6.0;
	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) {
			if (i == j) {
				K[N*i+j] = 0.0;                      // no "self-connections"!
			}
			else {
				K[N*i+j] = dt*(Kmean+Ksdev*randn()); // scale coupling constants by dt
			}
		}
	}

	// initialise oscillator phases with input (zero-mean Gaussian white noise)

	const double sqrtdt = sqrt(dt);
	const double Isdev = M_PI/40.0;  // noise intensity
	if (Isdev > 0.0) {
		for (size_t k=0; k<N*n; ++k) {
			h[k] = sqrtdt*Isdev*randn(); // scale input by sqrt(dt) [cf. Ornstein-Uhlenbeck process]
		}
	}
	else {
		// Note: if h allocated with calloc, it will be zeroed out;
		// if not, you should explicitly zero it out.
	}

	// integrate Kuramoto ODE

	const int RK4 = 0; // flag for RK4 (else Euler)
	if (RK4) {
		double* const kbuff = calloc(4*N,sizeof(double)); // see kuramoto_rk4()
		kuramoto_rk4(N,n,w,K,h,kbuff);
		free(kbuff);
	}
	else {
		kuramoto_euler(N,n,w,K,h);
	}

	// calculate order parameter

	order_param(N,n,h,r);

	// wrap oscillator phases to [-pi,pi) [if that's is what you want]
	//
	// phase_wrap(N*n,h);

	// write time stamp, order parameter and oscillator signals to file

	char ofile[] = "/tmp/kuramoto_demo.asc";     // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)(k+1)*dt);   // time stamp
		fprintf(fp," %17.8f",r[k]);              // order parameter
		for (size_t i=0; i<N; ++i) {
			fprintf(fp," %17.8f",sin(h[N*k+i])); // oscillator signal (waveform)
		}
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	// if Gnuplot installed (and accepts piped input - Linux, Mac),
	// display order parameter and oscillator signals. Else use your
	// favourite plotting program on data in output file.

#ifdef _GNUPLOT_HAVE_PIPE
	FILE* const gp = popen("gnuplot","w");
	if (gp == NULL) perror("failed to open pipe to Gnuplot\n");
	fprintf(gp,"set xlabel \"time\"\n");
	fprintf(gp,"set ylabel \"mean phase\"\n");
	fprintf(gp,"set key right bottom Left rev\n");
	fprintf(gp,"# set grid\n");
	fprintf(gp,"set xr [0:%g]\n",T);
	fprintf(gp,"set yr [0:1.05]\n");
	fprintf(gp,"set ytics 0.5\n");
	fprintf(gp,"set multiplot title \"Kuramoto demo\" layout 2,1\n");
	fprintf(gp,"set title \"Order parameter\"\n");
	fprintf(gp,"plot \"%s\" u 1:2 w l not\n",ofile);
	fprintf(gp,"set title \"Oscillator signals (waveforms)\"\n");
	fprintf(gp,"set ylabel \"amplitude\"\n");
	fprintf(gp,"set yr [-1.05:1.05]\n");
	fprintf(gp,"plot \\\n");
	for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,3+i);
	fprintf(gp,"NaN not\n");
	fprintf(gp,"unset multiplot\n");
	if (pclose(gp) == -1) perror("failed to close pipe to Gnuplot\n");
#endif

	// free memory

	free(r);
	free(h);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
