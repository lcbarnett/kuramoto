#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

// Program to demonstrate usage of Kuramoto C library.

int demo(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(N,      size_t, 4,         "number of oscillators");
	CLAP_CARG(T,      double, 1.0,       "total time (secs)");
	CLAP_CARG(fs,     double, 1000.0,    "Sampling frequency (Hz)");
	CLAP_CARG(wmax,   double, 20.0,      "oscillator frequency max (Hz)");
	CLAP_CARG(wmin,   double, 0.0,       "oscillator frequency min (Hz)");
	CLAP_CARG(Kmean,  double, 3.0,       "coupling constants mean (Hz)");
	CLAP_CARG(Ksdev,  double, Kmean/5.0, "coupling constants std. dev. (Hz)");
	CLAP_CARG(nmean,  double, 0.1,       "oscillator input noise magnitude mean (zero for no noise)");
	CLAP_CARG(nsdev,  double, nmean/5.0, "oscillator input noise magnitude std. dev.");
	CLAP_CARG(RK4,    int,    0,         "RK4 solver flag (else Euler)");
	CLAP_CARG(rseed,  ulong,  0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------");

	// seed random number generator (from command line if you want predictability)

	mt_t rng;
	mtuint_t seed = mt_seed(&rng,rseed);
	printf("\nrandom seed = %lu\n\n",seed);

	// some convenient constants

	const size_t n    = (size_t)ceil(T*fs); // number of integration steps
	const size_t m    = N*n;                // size of oscillator buffers
	const size_t M    = N*N;                // number of coupling constants
	const double ooN  = 1.0/(double)N;      // 1/N
	const double ffac = (2.0*M_PI)/fs;      // frequency scaling factor
	const double nfac = sqrt(ffac);         // Wiener noise scaling factor
	const double ffoN = ffac*ooN;           // coupling constants scaling factor

	// allocate memory

	double* const wdt = calloc(N,sizeof(double)); // oscillator frequencies
	double* const Kdt = calloc(M,sizeof(double)); // coupling constants
	double* const h   = calloc(m,sizeof(double)); // oscillator phases
	double* const r   = calloc(n,sizeof(double)); // order parameter
	double* const x   = calloc(m,sizeof(double)); // oscillator signal
	double* const y   = calloc(n,sizeof(double)); // oscillator agregated signal

	// random frequencies (uniform)

	for (size_t i=0; i<N; ++i) wdt[i] = ffac*(wmin+(wmax-wmin)*mt_rand(&rng));

	// random coupling constants (normal distribution), scaled by number of oscillators

	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) Kdt[N*i+j] = i == j ? 0.0 : ffoN*(Kmean+Ksdev*mt_randn(&rng)); // no "self-connections"!
	}

	// oscillator input noise (Wiener, with magnitudes log-normally distributed per oscillator)

	if (nmean > 0.0) {
		const double lnv = log(1.0+(nsdev*nsdev)/(nmean*nmean));
		const double mu  = log(nmean)-0.5*lnv;
		const double sig = sqrt(lnv);
		for (size_t i=0; i<N; ++i) {
			const double nmagi = nfac*exp(mu+sig*mt_randn(&rng)); // log-normal magnitude, Wiener scaling
			for (size_t k=i; k<m+i; k += N) h[k] = nmagi*mt_randn(&rng);
		}
	}
	else {
		memset(h,0,m*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	}

	// initial phases uniformly distributed on [0,2pi)

	for (size_t i=0; i<N; ++i) h[i] = 2.0*M_PI*mt_rand(&rng);


	// integrate Kuramoto ODE

	if (RK4) {
		double* const kbuff = calloc(4*N,sizeof(double)); // see kuramoto_rk4()
		kuramoto_rk4(N,n,wdt,Kdt,h,kbuff);
		free(kbuff);
	}
	else {
		kuramoto_euler(N,n,wdt,Kdt,h);
	}

	// calculate order parameter

	kuramoto_order_param(N,n,h,r,NULL);

	// wrap oscillator phases to [-pi,pi) [if that's what you want]
	//
	// phase_wrap(m,h);

	// generate signal from phases

	for (size_t j=0; j<m; ++j) x[j] = sin(h[j]);
	for (size_t k=0; k<n; ++k) {
		double yk = 0.0;
		for (size_t i=0; i<N; ++i) yk += x[N*k+i];
		y[k] = ooN*yk;
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
	free(Kdt);
	free(wdt);

	return EXIT_SUCCESS;
}
