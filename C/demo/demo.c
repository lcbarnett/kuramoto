#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "ode.h"

// Program to demonstrate usage of Kuramoto C library.

int demo(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(ode,    int,    1,         "1 - Euler, 2 - Heun, 3 - RK4");
	CLAP_CARG(N,      size_t, 4,         "number of oscillators");
	CLAP_CARG(T,      double, 1.0,       "total time (secs)");
	CLAP_CARG(fs,     double, 1000.0,    "Sampling frequency (Hz)");
	CLAP_CARG(wmax,   double, 20.0,      "oscillator frequency max (Hz)");
	CLAP_CARG(wmin,   double, 0.0,       "oscillator frequency min (Hz)");
	CLAP_CARG(Kmean,  double, 3.0,       "coupling constants mean (Hz)");
	CLAP_CARG(Ksdev,  double, Kmean/5.0, "coupling constants std. dev. (Hz)");
	CLAP_CARG(amean,  double, 0.0,       "phase lag mean (Hz)");
	CLAP_CARG(asdev,  double, amean/5.0, "phase lag mean std. dev. (Hz)");
	CLAP_CARG(nmean,  double, 0.1,       "oscillator input noise magnitude mean (sqrt(Hz) - zero for no noise)");
	CLAP_CARG(nsdev,  double, nmean/5.0, "oscillator input noise magnitude std. dev. (sqrt(Hz))");
	CLAP_CARG(RK4,    int,    0,         "RK4 solver flag (else Euler)");
	CLAP_CARG(rseed,  ulong,  0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,  "wxt",      "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------");

	const bool withpl = (amean > 0.0); // with phase lags
	const bool within = (nmean > 0.0); // with input noise

	// seed random number generator (from command line if you want predictability)

	mt_t rng;
	mtuint_t seed = mt_seed(&rng,rseed);
	printf("\nrandom seed = %lu\n\n",seed);

	// some convenient constants

	const size_t n    = (size_t)ceil(T*fs); // number of integration steps
	const size_t m    = N*n;                // size of oscillator buffers
	const size_t M    = N*N;                // number of coupling constants
	const double dt   = 1.0/fs;             // sample time step
	const double srdt = sqrt(dt);           // Wiener noise scaling factor

	// allocate memory

	double* const w = calloc(N,sizeof(double)); // oscillator frequencies
	double* const K = calloc(M,sizeof(double)); // coupling constants
	double* const a = withpl ? calloc(M,sizeof(double)) : NULL; // phase lags
	double* const h = calloc(m,sizeof(double)); // oscillator phases
	double* const r = calloc(n,sizeof(double)); // order parameter
	double* const x = calloc(m,sizeof(double)); // oscillator signal
	double* const y = calloc(n,sizeof(double)); // oscillator agregated signal

	// random frequencies (uniform)

	for (size_t i=0; i<N; ++i) w[i] = wmin+(wmax-wmin)*mt_rand(&rng);

	// random coupling constants (normal distribution), scaled by number of oscillators

	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) K[N*i+j] = i == j ? 0.0 : (Kmean+Ksdev*mt_randn(&rng))/(double)N; // no "self-connections"!
	}

	// phase lags (normal distribution)

	if (withpl) {
		for (size_t i=0; i<N; ++i) {
			for (size_t j=0; j<N; ++j) a[N*i+j] = amean + asdev*mt_randn(&rng);
		}
	}

	// oscillator input noise (Wiener, with magnitudes log-normally distributed per oscillator)

	if (within) {
		const double lnv = log(1.0+(nsdev*nsdev)/(nmean*nmean));
		const double mu  = log(nmean)-0.5*lnv;
		const double sig = sqrt(lnv);
		for (size_t i=0; i<N; ++i) {
			const double nmagi = exp(mu+sig*mt_randn(&rng)); // log-normal magnitude
			for (size_t k=i; k<m+i; k += N) h[k] = srdt*nmagi*mt_randn(&rng); // Wiener scaling
		}
	}
	else {
		memset(h,0,m*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	}

	// initial phases uniformly distributed on [-1,1)

	for (size_t i=0; i<N; ++i) h[i] = 2.0*mt_rand(&rng)-1.0;

	// integrate Kuramoto ODE

	const double ts = timer_start("Running ODE solver");
	if (withpl) {
		ODE(ode,kmotopl,h,N,n,dt,w,K,a);
	}
	else {
		ODE(ode,kmoto,h,N,n,dt,w,K);
	}
	timer_stop(ts);

	// calculate order parameter

	kuramoto_order_param(N,n,h,r,NULL);

	// generate signal from phases

	for (size_t j=0; j<m; ++j) x[j] = sin(TWOPI*h[j]);
	for (size_t k=0; k<n; ++k) {
		double yk = 0.0;
		for (size_t i=0; i<N; ++i) yk += x[N*k+i];
		y[k] = yk/(double)N;
	}

	// write time stamp, order parameter and oscillator signals to file

	char ofile[] = "/tmp/kuramoto_demo.asc";     // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)(k+1)/fs); // time stamp
		fprintf(fp," %17.8f",r[k]);            // order parameter
		fprintf(fp," %17.8f",y[k]);            // aggregate signal
		for (size_t i=0; i<N; ++i) fprintf(fp," %17.8f",x[N*k+i]); // signal
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
	fprintf(gp,"set term \"%s\" title \"Kuramoto oscillator demo\" size 1600,1200\n",gpterm);
	fprintf(gp,"datfile = \"%s\"\n",ofile);
	fprintf(gp,"set xlabel \"time\"\n");
	fprintf(gp,"set ylabel \"mean phase\"\n");
	fprintf(gp,"set key right bottom Left rev\n");
	fprintf(gp,"# set grid\n");
	fprintf(gp,"set xr [0:%g]\n",T);
	fprintf(gp,"set yr [0:1.05]\n");
	fprintf(gp,"set ytics 0.5\n");
	fprintf(gp,"set multiplot layout 3,1\n");
	fprintf(gp,"set title \"Order parameter\"\n");
	fprintf(gp,"plot datfile u 1:2 w l not\n");
	fprintf(gp,"set yr [-1.05:1.05]\n");
	fprintf(gp,"set title \"Aggregate oscillator signal (waveform)\"\n");
	fprintf(gp,"plot datfile u 1:3 w l not\n");
	fprintf(gp,"set title \"Oscillator signals (waveforms)\"\n");
	fprintf(gp,"set ylabel \"amplitude\"\n");
	fprintf(gp,"plot \\\n");
	for (size_t i=0; i<N; ++i) fprintf(gp,"datfile u 1:%zu w l not ,\\\n",4+i);
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

	free(y);
	free(x);
	free(r);
	free(h);
	if (withpl) free(a);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
