#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

// Program to demonstrate usage of Kuramoto C library.

int scratch(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(N,      size_t, 4,         "number of oscillators");
	CLAP_CARG(T,      double, 200.0,     "total integration time");
	CLAP_CARG(dt,     double, 0.01,      "integration step size");
	CLAP_CARG(wmean,  double, 0.0,       "oscillator frequencies mean (Hz)");
	CLAP_CARG(wsdev,  double, 0.2,       "oscillator frequencies std. dev. (Hz)");
	CLAP_CARG(Kmean,  double, 0.1,       "coupling constants mean (Hz)");
	CLAP_CARG(Ksdev,  double, Kmean/5.0, "coupling constants std. dev. (Hz)");
	CLAP_CARG(rseed,  ulong,  0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type (\"none\" for no plotting)");
#endif
	puts("---------------------------------------------------------------------------------------");

	// seed random number generator (from command line if you want predictability)

	mt_t rng;
	mtuint_t seed = mt_seed(&rng,rseed);
	printf("\nrandom seed = %lu\n\n",seed);

	// some convenient constants

	const size_t n   = (size_t)ceil(T/dt); // number of integration steps
	const size_t m   = N*n; // size of oscillator buffers
	const size_t M   = N*N; // number of coupling constants
	const double ooN = 1.0/(double)N;

	// allocate memory

	double* const wdt = calloc(N,sizeof(double)); // oscillator frequencies
	double* const Kdt = calloc(M,sizeof(double)); // coupling constants
	double* const h   = calloc(m,sizeof(double)); // oscillator phases
	double* const x   = calloc(m,sizeof(double)); // oscillator signal

	darray* const Kdtx = matalloc(N,N,NULL);
	darray* const hx   = matalloc(n,N,NULL);
	darray* const xx   = matalloc(n,N,NULL);

	// random frequencies (normal distribution)

	for (size_t i=0; i<N; ++i) wdt[i] = dt*TWOPI*(wmean+wsdev*mt_randn(&rng));

	// random coupling constants (normal distribution)

	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) {
			if (i == j) {
				Kdt[N*i+j] = 0.0; // no "self-connections"!
			}
			else {
				Kdt[N*i+j] = dt*TWOPI*ooN*(Kmean+Ksdev*mt_randn(&rng)); // scale coupling constants by N
			}
		}
	}

	memcpy(*Kdtx,Kdt,M*sizeof(double));

	// integrate Kuramoto ODE

	memset(h,0,    m*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	kuramoto_euler(N,n,wdt,Kdt,h); // dry run :-)

	memset(h,0,    m*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	const double ts1 = timer_start("algo 1");
	kuramoto_euler(N,n,wdt,Kdt,h);
	timer_stop(ts1);

	memset(*hx,0,m*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	const double ts2 = timer_start("algo 2");
	kuramoto_euler_alt(N,n,wdt,Kdtx,hx);
	timer_stop(ts2);

	// generate signal from phases

	for (size_t j=0; j<m; ++j) x[j] = sin(h[j]);

	for (size_t j=0; j<m; ++j) xx[0][j] = sin(hx[0][j]);


	// write time stamp, order parameter and oscillator signals to file

	char ofile[] = "/tmp/kuramoto_demo.asc";     // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t t=0; t<n; ++t) {
		fprintf(fp,"%17.8f",(double)(t+1)*dt); // time stamp
		for (size_t i=0; i<N; ++i) fprintf(fp," %17.8f",x[N*t+i]); // signal
		for (size_t i=0; i<N; ++i) fprintf(fp," %17.8f",xx[t][i]); // signal
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	if  (strncasecmp(gpterm,"none",4) != 0) {
		char gfile[] = "/tmp/kuramoto_demo.gp"; // Gnuplot command file
		FILE* const gp = fopen(gfile,"w");
		if (gp == NULL) {
			perror("failed to open Gnuplot command file\n");
			return EXIT_FAILURE;
		}
		fprintf(gp,"set term \"%s\" title \"Kuramoto oscillator demo\" size 1600,1200\n",gpterm);
		fprintf(gp,"set xlabel \"time\"\n");
		fprintf(gp,"set ylabel \"mean phase\"\n");
		fprintf(gp,"set key right bottom Left rev\n");
		fprintf(gp,"# set grid\n");
		fprintf(gp,"set xr [0:%g]\n",T);
		fprintf(gp,"set yr [-1.05:1.05]\n");
		fprintf(gp,"set ytics 0.5\n");
		fprintf(gp,"set ylabel \"amplitude\"\n");
		fprintf(gp,"set multiplot layout 2,1\n");
		fprintf(gp,"set title \"Oscillator signals\"\n");
		fprintf(gp,"plot \\\n");
		for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,i+2);
		fprintf(gp,"NaN not\n");
		fprintf(gp,"set title \"Oscillator signals x\"\n");
		fprintf(gp,"plot \\\n");
		for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,i+2+N);
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
	}
#endif

	// free memory

	matfree(xx);
	matfree(hx);
	matfree(Kdtx);

	free(x);
	free(h);
	free(Kdt);
	free(wdt);

	return EXIT_SUCCESS;
}
