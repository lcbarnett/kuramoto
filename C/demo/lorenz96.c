#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "ode.h"

// Program to demonstrate chaotic attractor system (Lorenz or Rossler).

int lorenz96(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type     default     description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(ode,    int,     1,         "1 - Euler, 2 - Heun, 3 - RK4");
	CLAP_CARG(N,      size_t,  7,         "number of variables");
	CLAP_VARG(T,      double,  1000.0,    "total integration time");
	CLAP_CARG(dt,     double,  0.01,      "integration step size");
	CLAP_VARG(F,      double,  8.0,       "Lorenz96 forcing parameter");
	CLAP_CARG(x0m,    double,  0.0,       "initial value for variables - mean");
	CLAP_CARG(x0s,    double,  2.0,       "initial value for variables - std. dev.");
	CLAP_CARG(nmag,   double,  0.0,       "noise magnitude");
	CLAP_CARG(rseed,  uint,    0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,   "wxt",      "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------");

	if (N < 4) {
		fprintf(stderr,"\nMust be at least 4 variables\n\n");
		return EXIT_FAILURE;
	}

	// which ODE solver?

	if (ode < 1 || ode > 3) {
		fprintf(stderr,"\nUnknown ODE solver\n\n");
		return EXIT_FAILURE;
	}

	// seed random number generator (from command line if you want predictability)

	mt_t rng;
	mtuint_t seed = mt_seed(&rng,rseed);
	printf("\nrandom seed = %lu\n\n",seed);

	// some convenient constants

	const size_t n = (size_t)ceil(T/dt); // number of integration steps
	T = dt*(double)n;
	printf("time increments = %lu\n\n",n);

	// allocate memory

	const size_t m = N*n;
	double* const x = calloc(m,sizeof(double)); // 3D variable

	// noise

	if (nmag > 0.0) {
		for (size_t i=0; i<m; ++i) x[i] = sqrt(dt)*nmag*mt_randn(&rng);
	}
	else {
		memset(x,0,m*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
	}

	// initialise

	for (size_t k=0; k<N; ++k) x[k] = x0m + x0s*mt_randn(&rng);

	// integrate Lorenz96 ODE

	const double ts = timer_start("Running ODE solver");
	ODE(ode,lrnz96,x,N,n,dt,F);
	timer_stop(ts);

	// write variables to file

	char ofile[] = "/tmp/chaos_demo.asc"; // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; k += N) {
		fprintf(fp,"%17.8f %17.8f %17.8f\n",x[k],x[k+1],x[k+2]);
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	// finished with x

	free(x);

	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	char gfile[] = "/tmp/chaos_demo.gp"; // Gnuplot command file
	FILE* const gp = fopen(gfile,"w");
	if (gp == NULL) {
		perror("failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gp,"set term \"%s\" size 1600,1200\n",gpterm);
	fprintf(gp,"set mouse ruler\n");
	fprintf(gp,"unset key\n");
	fprintf(gp,"set grid\n");
	fprintf(gp,"set title \"Lorenz 96 system (%s) : %zu variables, F = %g\"\n",ode == 1 ? "Euler" : ode == 2 ? "Heun" : ode == 3 ? "RK4" : "", N,F);
	fprintf(gp,"set xlabel \"x\"\n");
	fprintf(gp,"set ylabel \"y\"\n");
	fprintf(gp,"set zlabel \"z\"\n");
	fprintf(gp,"splot \"%s\" u 1:2:3 w l not\n",ofile);
	if (fclose(gp) != 0) {
		perror("Failed to close Gnuplot command file");
		return EXIT_FAILURE;
	}
	const size_t strlen = 100;
	char gpcmd[strlen+1];
	snprintf(gpcmd,strlen,"gnuplot -p %s",gfile);
	printf("Gnuplot command: %s\n\n",gpcmd);
	if (system(gpcmd) == -1) {
		perror("Failed to run Gnuplot command");
		return EXIT_FAILURE;
	}
#endif

	return EXIT_SUCCESS;
}
