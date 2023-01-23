#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "kuramoto.h"

// Uniform random double on [0,1) [Note: you might want a better PRNG]

static inline double randu()
{
	return (double)random()/((double)(RAND_MAX)+1.0);
}

// Standard normal random double (Box-Muller)

static inline double randn()
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

// Main function. Call as:
//
// kuramoto_demo <seed>

int main(int argc, char *argv[])
{
	// seed random number generator

	srand(argc > 1 ? (unsigned)atoi(argv[1]) : 1);

	// Kuramoto model size parameters

	const size_t N  = 4;                          // number of oscillators
	const double T  = 100.0;                      // total integration time
	const double dt = 0.01;                       // integration step size
	const size_t n  = (size_t)round(T/dt);        // number of integration steps

	// allocate memory

	double* const w = calloc(N,  sizeof(double)); // oscillator frequencies
	double* const K = calloc(N*N,sizeof(double)); // coupling constants
	double* const h = calloc(N*n,sizeof(double)); // oscillator phases, unwrapped
	double* const p = calloc(N*n,sizeof(double)); // oscillator phases. wrapped to [-pi,pi)
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
	const double Isdev = M_PI/20.0;
	for (size_t k=0; k<N*n; ++k) {
		h[k] = sqrtdt*Isdev*randn(); // scale input by sqrt(dt) [cf. Ornstein-Uhlenbeck process]
	}

	// integrate Kuramoto ODE

	kuramoto_euler(N,n,w,K,h);

	// calculate order parameter

	order_param(N,n,h,r);

	// wrap oscillator phases to [-pi,pi) [if this is what you want]

	phase_wrap(N,n,h,p);

	// write time stamp, order parameter and oscillator phases to file

	FILE* const fp = fopen("/tmp/kuramoto_demo.asc","w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)(k+1)*dt); // time stamp
		fprintf(fp," %17.8f",r[k]);            // order parameter
		for (size_t i=0; i<N; ++i) {
			fprintf(fp," %17.8f",p[N*k+i]);   // oscillator phase (wrapped)
		}
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	// free memory

	free(r);
	free(p);
	free(h);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
