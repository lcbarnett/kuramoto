#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

// Program to demonstrate chaotic attractor system (Lorenz or Rossler).

int chaos(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type     default     description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_VARG(T,      double,  500.0,     "total integration time");
	CLAP_CARG(dt,     double,  0.01,      "integration step size");
	CLAP_CARG(csys,   cstr,   "lorenz",   "chaotic system (lorenz or rossler)");
//	CLAP_CARG(RK4,    int,     0,         "RK4 solver flag (else Euler)");
	CLAP_VARG(p1,     double,  NAN,       "parameter 1 (NAN for system default)");
	CLAP_VARG(p2,     double,  NAN,       "parameter 2 (NAN for system default)");
	CLAP_VARG(p3,     double,  NAN,       "parameter 3 (NAN for system default)");
	CLAP_CARG(nmag,   double,  0.0,       "noise magnitude");
	CLAP_CARG(rseed,  uint,    0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,    GPTERM,    "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------");

	// seed random number generator (from command line if you want predictability)

	mt_t rng;
	mtuint_t seed = mt_seed(&rng,rseed);
	printf("\nrandom seed = %lu\n\n",seed);

	// some convenient constants

	const size_t n = (size_t)ceil(T/dt); // number of integration steps
	T = dt*(double)n;
	printf("time increments = %lu\n\n",n);

	// allocate memory

	const size_t N = 3*n;
	double* const x = calloc(N,sizeof(double)); // 3D variable

	// noise

	if (nmag > 0.0) {
		for (size_t i=0; i<N; ++i) x[i] = sqrt(dt)*nmag*mt_randn(&rng);
	}
	else {
		memset(x,0,N*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
	}

	if (strncasecmp(csys,"lorenz",6) == 0) {
		if (isnan(p1)) p1 = 10.0;
		if (isnan(p2)) p2 = 28.0;
		if (isnan(p3)) p3 = 8.0/3.0;
		lorenz_euler(n,dt,p1,p2,p3,x);
	}
	else if (strncasecmp(csys,"rossler",7) == 0) {
		if (isnan(p1)) p1 = 0.1;
		if (isnan(p2)) p2 = 0.1;
		if (isnan(p3)) p3 = 14.0;
		rossler_euler(n,dt,p1,p2,p3,x);
	}
	else {
		fprintf(stderr,"Unknown chaotic system \"%s\"",csys);
		return EXIT_FAILURE;
	}

	// write time stamp and variables to file

	char ofile[] = "/tmp/chaos_demo.asc";     // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		const size_t r = 3*k;
		fprintf(fp,"%17.8f %17.8f %17.8f %17.8f\n",dt*(double)(k+1),x[r],x[r+1],x[r+2]);            // order parameter
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}
/*
	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	char gfile[] = "/tmp/stulan_demo.gp"; // Gnuplot command file
	FILE* const gp = fopen(gfile,"w");
	if (gp == NULL) {
		perror("failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gp,"set term \"%s\" title \"Coupled Stuart-Landau oscillators\" size 1600,1200\n",gpterm);
	fprintf(gp,"set xlabel \"time\"\n");
	fprintf(gp,"set key right bottom Left rev\n");
	fprintf(gp,"# set grid\n");
//	fprintf(gp,"set xr [0:%g]\n",T);
	fprintf(gp,"set yr [0:*]\n");
//	fprintf(gp,"set ytics 0.5\n");
	fprintf(gp,"set multiplot layout 3,1\n");
	fprintf(gp,"set title \"Order parameter\"\n");
	fprintf(gp,"set ylabel \"mean magnitude\"\n");
	fprintf(gp,"plot \"%s\" u 1:2 w l not\n",ofile);
	fprintf(gp,"set title \"Oscillator magnitudes\"\n");
	fprintf(gp,"set ylabel \"magnitude\"\n");
	fprintf(gp,"plot \\\n");
	for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,3+i);
	fprintf(gp,"NaN not\n");
	fprintf(gp,"set title \"Oscillator signals\"\n");
	fprintf(gp,"set ylabel \"amplitude\"\n");
	fprintf(gp,"set yr [*:*]\n");
	fprintf(gp,"plot \\\n");
	for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" u 1:%zu w l not ,\\\n",ofile,3+N+i);
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
*/
	free(x);

	return EXIT_SUCCESS;
}
