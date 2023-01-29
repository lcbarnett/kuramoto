#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "kuramoto.h"

// Program to demonstrate usage of Kuramoto C library.

int demo(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches; e.g.:
	//
	// kuramoto_demo -N 10 -T 1000 -dt 0.001 -Isdev 0
	//
	// Arg:  name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_ARG(N,      size_t, 4,         "number of oscillators");
	CLAP_ARG(T,      double, 200.0,     "total integration time");
	CLAP_ARG(dt,     double, 0.01,      "integration step size");
	CLAP_ARG(wmean,  double, 0.0,       "oscillator frequencies mean");
	CLAP_ARG(wsdev,  double, M_PI/7.0,  "oscillator frequencies std. dev.");
	CLAP_ARG(Kmean,  double, 0.8/N,     "coupling constants mean");
	CLAP_ARG(Ksdev,  double, Kmean/6.0, "coupling constants std. dev.");
	CLAP_ARG(Isdev,  double, M_PI/40.0, "input noise intensity (zero for deterministic)");
	CLAP_ARG(rseed,  uint,   0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_ARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type (if available)");
#endif
	puts("---------------------------------------------------------------------------------------");

	// seed random number generator (from command line if you want predictability)

	const uint seed = rseed > 0 ? rseed : get_rand_seed();
	printf("\nrandom seed = %u\n\n",seed);
	srand(seed);

	// number of integration steps

	const size_t n  = (size_t)ceil(T/dt);         // number of integration steps

	// allocate memory

	double* const w = calloc(N,  sizeof(double)); // oscillator frequencies
	double* const K = calloc(N*N,sizeof(double)); // coupling constants
	double* const h = calloc(N*n,sizeof(double)); // oscillator phases, unwrapped
	double* const r = calloc(n,  sizeof(double)); // order parameter

	// random frequencies (normal distribution)

	for (size_t i=0; i<N; ++i) {
		w[i] = dt*(wmean+wsdev*randn()); // scale frequencies by dt
	}

	// random coupling constants (normal distribution)

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
	if (Isdev > 0.0) {
		for (size_t k=0; k<N*n; ++k) {
			h[k] = sqrtdt*Isdev*randn(); // scale input by sqrt(dt) [cf. Ornstein-Uhlenbeck process]
		}
	}
	else {
		memset(h,0,N*N*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
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

	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	char gfile[] = "/tmp/kuramoto_demo.gp"; // Gnuplot command file
	FILE* const gp = fopen(gfile,"w");
	if (gp == NULL) {
		perror("failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gp,"set term \"%s\" title \"Kuramoto oscillator demo\" size 1600,800\n",gpterm);
	fprintf(gp,"set xlabel \"time\"\n");
	fprintf(gp,"set ylabel \"mean phase\"\n");
	fprintf(gp,"set key right bottom Left rev\n");
	fprintf(gp,"# set grid\n");
	fprintf(gp,"set xr [0:%g]\n",T);
	fprintf(gp,"set yr [0:1.05]\n");
	fprintf(gp,"set ytics 0.5\n");
	fprintf(gp,"set multiplot layout 2,1\n");
	fprintf(gp,"set title \"Order parameter\"\n");
	fprintf(gp,"plot \"%s\" u 1:2 w l not\n",ofile);
	fprintf(gp,"set title \"Oscillator signals (waveforms)\"\n");
	fprintf(gp,"set ylabel \"amplitude\"\n");
	fprintf(gp,"set yr [-1.05:1.05]\n");
	fprintf(gp,"plot \\\n");
	for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,3+i);
	fprintf(gp,"NaN not\n");
	fprintf(gp,"unset multiplot\n");
	if (fclose(gp) != 0) {
		perror("Failed to close Gnuplot command file");
		return EXIT_FAILURE;
	}
	const size_t strlen = 100;
	char gpcmd[strlen+1];
	strncpy(gpcmd,"gnuplot -p ",strlen);
	strncat(gpcmd,gfile,strlen);
	printf("Gnuplot command: %s\n\n",gpcmd);
	if (system(gpcmd) == -1) {
		perror("Failed to run Gnuplot command");
		return EXIT_FAILURE;
	}
#endif

	// free memory

	free(r);
	free(h);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
