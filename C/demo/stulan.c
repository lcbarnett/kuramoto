#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "mt64.h"
#include "kuramoto.h"

// Program to demonstrate usage of coupled Stuart-Landau oscillators.

int stulan(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches; e.g.:
	//
	// kuramoto stulan -N 10 -T 1000 -dt 0.001 -Isdev r0 = 3
	//
	// Arg:   name    type    default    description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(N,      size_t, 6,         "number of oscillators");
	CLAP_CARG(T,      double, 500.0,     "total integration time");
	CLAP_CARG(dt,     double, 0.01,      "integration step size");
	CLAP_CARG(wmean,  double, 0.0,       "oscillator frequencies mean");
	CLAP_CARG(wsdev,  double, 0.3,       "oscillator frequencies std. dev.");
	CLAP_CARG(Kmean,  double, 0.1,       "coupling constants mean");
	CLAP_CARG(Ksdev,  double, Kmean/8.0, "coupling constants std. dev.");
	CLAP_CARG(amean,  double, 1.0,       "growth constants mean");
	CLAP_CARG(asdev,  double, amean/8.0, "growth constants std. dev.");
	CLAP_CARG(Isdev,  double, 0.0,       "input noise intensity");
//	CLAP_CARG(RK4,    int,    0,         "RK4 solver flag (else Euler)");
	CLAP_CARG(rseed,  uint,   0,         "random seed (or 0 for random random seed)");
#ifdef _HAVE_GNUPLOT
	CLAP_CARG(gpterm, cstr,   GPTERM,    "Gnuplot terminal type (if available)");
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

	double* const w = calloc(N,sizeof(double)); // oscillator frequencies
	double* const K = calloc(M,sizeof(double)); // oscillator coupling constants
	double* const a = calloc(N,sizeof(double)); // oscillator growth constants
	double* const x = calloc(m,sizeof(double)); // oscillator real part
	double* const y = calloc(m,sizeof(double)); // oscillator imag part
	double* const r = calloc(n,sizeof(double)); // order parameter
	double* const s = calloc(m,sizeof(double)); // oscillator magnitudes

	// random frequencies (normal distribution)

	for (size_t i=0; i<N; ++i) w[i] = wmean+wsdev*mt_randn(&rng);

	// random coupling constants (normal distribution)

	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) {
			if (i == j) {
				K[N*i+j] = 0.0; // no "self-connections"!
			}
			else {
				K[N*i+j] = ooN*(Kmean+Ksdev*mt_randn(&rng)); // scale coupling constants by N
			}
		}
	}

	// random growth constants (log-normal distribution)

	const double lna = log(amean*amean+asdev*asdev);
	const double lnb = 2.0*log(amean);
	const double am  = lnb-lna/2.0;
	const double as  = sqrt(lna-lnb);
	for (size_t i=0; i<N; ++i)  a[i] = exp(am+as*mt_randn(&rng));

	// input zero-mean Gaussian white noise

	if (Isdev > 0.0) {
		for (size_t k=0; k<m; ++k) x[k] = Isdev*mt_randn(&rng);
		for (size_t k=0; k<m; ++k) y[k] = Isdev*mt_randn(&rng);
	}
	else {
		memset(x,0,m*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
		memset(y,0,m*sizeof(double));  // zero-fill for no input [in fact here calloc will have done that]
	}

	// initialise oscillators (uniform random on unit circle)

//	const double r0adj = r0/(1.0+SLMAGIC*(double)(N-1));
	for (size_t i=0; i<N; ++i) {
		const double h = TWOPI*mt_rand(&rng);
		x[i] = cos(h);
		y[i] = sin(h);
	}

	// integrate Kuramoto ODE

//	if (RK4) {
//		double* const kbuff = calloc(4*N,sizeof(double)); // see kuramoto_rk4()
//		kuramoto_rk4(N,n,w,K,h,kbuff);
//		free(kbuff);
//	}
//	else {
		stulan_euler(N,n,dt,w,K,a,x,y);
//	}

	// calculate order parameter

	stulan_order_param(N,n,x,y,r,NULL);

	// magnitudes

	stulan_magnitudes(N,n,x,y,s);

	// write time stamp, order parameter and oscillator signals to file

	char ofile[] = "/tmp/stulan_demo.asc";     // output file (ASCII)
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)(k+1)*dt); // time stamp
		fprintf(fp," %17.8f",r[k]);            // order parameter
		for (size_t i=0; i<N; ++i) fprintf(fp," %17.8f",s[N*k+i]); // oscillator magnitudes
		for (size_t i=0; i<N; ++i) fprintf(fp," %17.8f",x[N*k+i]); // oscillator signals
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

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

	// free memory

	free(s);
	free(r);
	free(y);
	free(x);
	free(a);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
