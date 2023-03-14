#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "kuramoto.h"

// Program to demonstrate usage of coupled Stuart-Landau oscillators.

int scratch(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:  name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_ARG(nmax,   size_t, 40,        "maximum number of oscillators");
	CLAP_ARG(S,      size_t, 10000,     "number of samples");
	CLAP_ARG(rseed,  uint,   0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_ARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type (if available)");
#endif
	puts("---------------------------------------------------------------------------------------");

	const uint seed = rseed > 0 ? rseed : get_rand_seed();
	printf("\nrandom seed = %u\n\n",seed);
	srand(seed);

	double* const z = calloc(nmax,sizeof(double));

	for (size_t n=1; n<=nmax; ++n) {
		printf("n = %3zu : ",n);
		double zn = 0.0;
		for (size_t s=0; s<S; ++s) {
			double x = 0.0, y = 0.0;
			for (size_t i=0; i<n; ++i) {
				const double h = TWOPI*randu();
				x += cos(h);
				y += sin(h);
			}
			zn += hypot(x,y);
		}
		z[n-1] = zn/(double)S;
		printf("z = %6.2f\n",z[n-1]);
	}
	z[0] = 1.0; // because it really is :-)

	double fac1 = 0.0;
	for (size_t n=2; n<=nmax; ++n) fac1 += ((double)(n-1))*(z[n-1]-1.0);
	double fac2 = 0.0;
	for (size_t n=2; n<=nmax; ++n) fac2 += ((double)(n-1))*((double)(n-1));
	const double slope = fac1/fac2;

	printf("\nslope = %.12f\n",slope);

	char ofile[] = "/tmp/stulan_scratch.asc";
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t n=1; n<=nmax; ++n) {
		fprintf(fp,"%6zu",n);
		fprintf(fp," %17.8f",z[n-1]);
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	char gfile[] = "/tmp/stulan_scratch.gp"; // Gnuplot command file
	FILE* const gp = fopen(gfile,"w");
	if (gp == NULL) {
		perror("failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gp,"set term \"%s\" size 1600,1200\n",gpterm);
	fprintf(gp,"set title \"Oscillator mean amplitude scaling\"\n");
	fprintf(gp,"set xlabel \"number of oscillators\"\n");
	fprintf(gp,"set ylabel \"mean amplitude\"\n");
	fprintf(gp,"set key right bottom Left rev\n");
	fprintf(gp,"set grid\n");
	fprintf(gp,"set xr [0.5:%g]\n",(double)nmax+0.5);
	fprintf(gp,"set arrow from first 1,first 1 to first %zu,first %g nohead\n",nmax,1.0+((double)(nmax-1))*slope);
	fprintf(gp,"plot \"%s\" u 1:2 w p pt 7 not\n",ofile);
	if (fclose(gp) != 0) {
		perror("Failed to close Gnuplot command file");
		return EXIT_FAILURE;
	}
	const size_t strlen = 100;
	char gpcmd[strlen+1];
	strncpy(gpcmd,"gnuplot -p ",strlen);
	strncat(gpcmd,gfile,strlen);
	printf("\nGnuplot command: %s\n\n",gpcmd);
	if (system(gpcmd) == -1) {
		perror("Failed to run Gnuplot command");
		return EXIT_FAILURE;
	}
#endif

	// free memory

	free(z);

	return EXIT_SUCCESS;
}
