#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

typedef struct {
	size_t tnum;
	size_t nipert;
	size_t S;
	double* z;
} targ_t;

void* compfun(void *arg)
{
	targ_t* targ = (targ_t*)arg;

	mt_t rng;
	mt_seed(&rng,0);
	printf("thread %zu start\n",targ->tnum);
	fflush(stdout);

	const size_t toff = targ->tnum*targ->nipert+1;

	for (size_t n=toff; n<toff+targ->nipert; ++n) {
		const double sqrtn = sqrt((double)n);
		double zn = 0.0;
		for (size_t s=0; s<targ->S; ++s) {
			double x = 0.0, y = 0.0;
			for (size_t i=0; i<n; ++i) {
				const double theta = TWOPI*mt_rand(&rng);
				x += cos(theta);
				y += sin(theta);
			}
			zn += hypot(x,y);
		}
		targ->z[n] += zn/(sqrtn*(double)targ->S);
		printf("thread %zu : n = %3zu : z = %.12f\n",targ->tnum,n,targ->z[n]);
		fflush(stdout);
	}

	printf("thread %zu exit\n",targ->tnum);
	fflush(stdout);
	pthread_exit(NULL);
}

int scratch(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches.
	//
	// Arg:   name      type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(nmax,     size_t, 60,        "maximum number of oscillators");
	CLAP_CARG(nthreads, size_t, 3,         "maximum number of oscillators");
	CLAP_CARG(S,        size_t, 1000000,   "number of samples");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm,   cstr,   GPTERM,    "Gnuplot terminal type");
#endif
	puts("---------------------------------------------------------------------------------------\n");

	const size_t nipert = nmax/nthreads; // number of iterations per thread
	if (nthreads*nipert != nmax) {
		fprintf(stderr,"Error: number of threads must divide maximum number of oscillators\n");
		exit(EXIT_FAILURE);
	}

	double* const z = calloc(nmax+1,sizeof(double));

	pthread_t threads[nthreads];
	pthread_attr_t attr;

	// initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	// kick off computation threads
	targ_t targ[nipert];
	for (size_t tnum=0; tnum<nthreads; ++tnum) {
		targ[tnum].tnum = tnum;
		targ[tnum].nipert = nipert;
		targ[tnum].S = S;
		targ[tnum].z = z;
		const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targ[tnum]);
		if (tres) {
			fprintf(stderr,"Error: unable to create thread %zu",tnum);
			exit(EXIT_FAILURE);
		}
	}

	// free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for (size_t tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_join(threads[tnum],NULL);
		if (tres) {
			fprintf(stderr,"Error: unable to join thread %zu",tnum);
			exit(EXIT_FAILURE);
		}
	}

	z[1] = 1.0; // because it is :-)

	char ofile[] = "/tmp/stulan_scratch.asc";
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t n=1; n<=nmax; ++n) {
		fprintf(fp,"%6zu",n);
		fprintf(fp," %17.8f",z[n]);
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
	fprintf(gp,"set yr [0:*]\n");
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

	pthread_exit(NULL);
}
