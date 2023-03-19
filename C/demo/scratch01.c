#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

// Program to demonstrate usage of coupled Stuart-Landau oscillators.

int scratch(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:  name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_ARG(S,      size_t, 10000000,  "number of samples");
	CLAP_ARG(B,      size_t, 100,       "number of bins");
	CLAP_ARG(rseed,  ulong,  0,         "random seed (or 0 for random random seed)");
	CLAP_ARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type (if available)");
	puts("---------------------------------------------------------------------------------------");

	mt_t rng;
	mtuint_t seed1 = mt_seed(&rng,rseed);

	const uint seed2 = rseed > 0 ? rseed : get_rand_seed();
	srand(seed2);

	printf("\nrandom seed 1 = %zu\n",seed1);
	printf("random seed 2 = %u\n",seed2);

	const double DB = (double)B;
	const double DS = (double)S;

	double* const c1 = calloc(B,sizeof(size_t));
	double* const c2 = calloc(B,sizeof(size_t));

	for (size_t s=0; s<S; ++s) ++c1[(size_t)floor(DB*mt_rand(&rng))];
	for (size_t s=0; s<S; ++s) ++c2[(size_t)floor(DB*randu())];

	char ofile[] = "/tmp/stulan_scratch.asc";
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t b=0; b<B; ++b) {
		fprintf(fp,"%6zu",b);
		fprintf(fp," %17.8f",(double)(c1[b])/DS);
		fprintf(fp," %17.8f",(double)(c2[b])/DS);
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
	fprintf(gp,"set title \"Uniform rand histogram\"\n");
	fprintf(gp,"unset key\n");
	fprintf(gp,"set grid\n");
//	fprintf(gp,"set xr [0.5:%g]\n",(double)nmax+0.5);
	fprintf(gp,"set yr [0:*]\n");
	fprintf(gp,"plot \"%s\" u 1:2 w p pt 7 not, \\\n",ofile);
	fprintf(gp,"\"%s\" u 1:3 w p pt 7 not\n",ofile);
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

	free(c2);
	free(c1);

	return EXIT_SUCCESS;
}
