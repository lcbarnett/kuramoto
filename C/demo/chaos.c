#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "ode.h"

// Program to demonstrate chaotic attractor system (Lorenz or Rossler).

int chaos(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name    type     default     description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_VARG(T,      double,  200.0,     "total integration time");
	CLAP_CARG(dt,     double,  0.01,      "integration step size");
	CLAP_CARG(sys,    int,     LORENZ,    "system type : 0 - LORENZ, 1 - ROSSLER, 2 - THOMAS");
	CLAP_CARG(ode,    int,     EULER,     "ode solver  : 0 - EULER , 1 - HEUN,    2 - RKFOUR");
	CLAP_VARG(p1,     double,  NAN,       "parameter 1     (NAN for system default)");
	CLAP_VARG(p2,     double,  NAN,       "parameter 2     (NAN for system default)");
	CLAP_VARG(p3,     double,  NAN,       "parameter 3     (NAN for system default)");
	CLAP_VARG(x1,     double,  NAN,       "initial value 1 (NAN for system default)");
	CLAP_VARG(x2,     double,  NAN,       "initial value 2 (NAN for system default)");
	CLAP_VARG(x3,     double,  NAN,       "initial value 3 (NAN for system default)");
	CLAP_CARG(nmag,   double,  0.0,       "noise magnitude");
	CLAP_CARG(rseed,  uint,    0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,   "wxt",      "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------");

	const char* odename = ode2str(ode);
	if (odename == NULL) error(EXIT_FAILURE,0,"Unknown ODE solver");

	const char* sysname = sys2str(sys);
	if (sysname == NULL) error(EXIT_FAILURE,0,"Unknown system type");
	if (sys > THOMAS)    error(EXIT_FAILURE,0,"Unhandled system type");

	// parameter defaults

	switch (sys) {
		case LORENZ:
			if (isnan(p1)) p1 = 10.0;
			if (isnan(p2)) p2 = 28.0;
			if (isnan(p3)) p3 = 8.0/3.0;
			if (isnan(x1)) x1 = 1.0;
			if (isnan(x2)) x2 = 1.0;
			if (isnan(x3)) x3 = 1.0;
			break;
		case ROSSLER:
			if (isnan(p1)) p1 = 0.1;
			if (isnan(p2)) p2 = 0.1;
			if (isnan(p3)) p3 = 14.0;
			if (isnan(x1)) x1 = 1.0;
			if (isnan(x2)) x2 = 1.0;
			if (isnan(x3)) x3 = 1.0;
			break;
		case THOMAS:
			if (isnan(p1)) p1 = 0.1;
			if (isnan(x1)) x1 = 1.0;
			if (isnan(x2)) x2 = 2.0;
			if (isnan(x3)) x3 = 3.0;
			break;
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

	const size_t N = 3*n;
	double* const x = calloc(N,sizeof(double)); // 3D variable

	// noise

	if (nmag > 0.0) {
		for (size_t i=0; i<N; ++i) x[i] = sqrt(dt)*nmag*mt_randn(&rng);
	}
	else {
		memset(x,0,N*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
	}

	// initialise

	x[0] = x1;
	x[1] = x2;
	x[2] = x3;

	// run

	const double ts = timer_start("Running ODE solver");
	switch (sys) {
		case LORENZ  : ODE(ode,lorenz, x,3,n,dt,p1,p2,p3); break;
		case ROSSLER : ODE(ode,rossler,x,3,n,dt,p1,p2,p3); break;
		case THOMAS  : ODE(ode,thomas, x,3,n,dt,p1);       break;
	}
	timer_stop(ts);

	// write variables to file

	char ofile[] = "/tmp/chaos_demo.asc"; // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		const size_t r = 3*k;
		fprintf(fp,"%17.8f %17.8f %17.8f\n",x[r],x[r+1],x[r+2]);
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
	fprintf(gp,"set title \"%s system (%s)\"\n",sysname,odename);
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
