#include <stdio.h>
#include <string.h>
#include <math.h>

#include "clap.h"
#include "kutils.h"
#include "kuramoto.h"

#define CD_SRATE 44100.0
#define MIDDLE_C 261.625565

// Program to demonstrate usage of Kuramoto C library - as an audio synth :-)

int audio(int argc, char *argv[])
{
	// CLAP (command-line argument parser). Default values may
	// be overriden on the command line as switches; e.g.:
	//
	// kuramoto audio -N 10 -T 5 -f 20000 -Isdev 0
	//
	// Arg:  name    type    default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_ARG(N,      size_t, 6,            "number of oscillators");
	CLAP_ARG(T,      double, 5.0,          "total time (seconds)");
	CLAP_ARG(f,      double, CD_SRATE,     "sampling frequency (Hz)");
	CLAP_ARG(wmean,  double, 0.0,          "oscillator frequencies mean (Hz)");
	CLAP_ARG(wsdev,  double, 2.0*MIDDLE_C, "oscillator frequencies std. dev. (Hz)");
	CLAP_ARG(Kmean,  double, 4.0,          "coupling constants mean (dimensionless)");
	CLAP_ARG(Ksdev,  double, Kmean/6.0,    "coupling constants std. dev. (dimensionless)");
	CLAP_ARG(Kbias,  double, 0.5,          "coupling constants bias (probability");
	CLAP_ARG(Isdev,  double, 0.2,          "input noise intensity (Hz: zero for deterministic)");
	CLAP_ARG(RK4,    int,    0,            "RK4 solver flag (else Euler)");
	CLAP_ARG(rseed,  uint,   0,            "random seed (or 0 for random random seed)");
	CLAP_ARG(pcm,    int,    16,           "PCM bits: 16 or 24 (unsigned), -32 or -64 (fp), or zero for no PCM");
	CLAP_ARG(cagg,   int,    1,            "PCM aggregate channels?");
#ifdef _HAVE_GNUPLOT
	CLAP_ARG(Ts,     double, 1.0,          "display time start (seconds)");
	CLAP_ARG(Te,     double, 1.1,          "display time end   (seconds)");
	CLAP_ARG(gpterm, cstr,   GPTERM,       "Gnuplot terminal type (if available) or \"noplot\"");
#endif
	puts("---------------------------------------------------------------------------------------");

	// seed random number generator (from command line if you want predictability)

	const uint seed = rseed > 0 ? rseed : get_rand_seed();
	printf("\nrandom seed = %u\n\n",seed);
	srand(seed);

	// some convenient constants

	const double dt  = 1.0/f;              // integration time step size (1/sampling frequency should do)
	const int    F   = (int)round(f);      // nearest int to sampling frequency
	const size_t n   = (size_t)ceil(T/dt); // number of integration steps
	const size_t m   = N*n;                // size of oscillator buffers
	const size_t M   = N*N;                // number of coupling constants
	const double ooN = 1.0/(double)N;

	// allocate memory

	double* const w = calloc(N,sizeof(double)); // oscillator frequencies
	double* const K = calloc(M,sizeof(double)); // coupling constants
	double* const h = calloc(m,sizeof(double)); // oscillator phases
	double* const r = calloc(n,sizeof(double)); // order parameter
	double* const x = calloc(m,sizeof(double)); // oscillator signal
	double* const y = calloc(n,sizeof(double)); // oscillator agregated signal

	// random frequencies (normal distribution)

	for (size_t i=0; i<N; ++i) {
		w[i] = TWOPI*(wmean+wsdev*randn());
	}

	// random coupling constants (normal distribution); note that we take incoming couplings as
	// multipliers for the corresponding oscillator's natural frequency - so that the become
	// dimensionless - and scale by the number of oscillators

	for (size_t i=0; i<N; ++i) {
		const double ooNwi = ooN*w[i]; // multiplier is w[i]/N
		for (size_t j=0; j<N; ++j) {
			if (i == j) {
				K[N*i+j] = 0.0; // no "self-connections"!
			}
			else {
				K[N*i+j] = ooNwi*((randu()<Kbias?Kmean:-Kmean)+Ksdev*randn()); // scale coupling constants by w[i] and N
			}
		}
	}

	// initialise oscillator phases with input (zero-mean Gaussian white noise)

	if (Isdev > 0.0) {
		for (size_t k=0; k<m; ++k) {
			h[k] = TWOPI*Isdev*randn();
		}
	}
	else {
		memset(h,0,M*sizeof(double)); // zero-fill for no input [in fact here calloc will have done that]
	}

	// integrate Kuramoto ODE

	const double ts1 = timer_start("simulating Kuramoto system");
	if (RK4) {
		double* const kbuff = calloc(4*N,sizeof(double)); // see kuramoto_rk4()
		kuramoto_rk4(N,n,dt,w,K,h,kbuff);
		free(kbuff);
	}
	else {
		kuramoto_euler(N,n,dt,w,K,h);
	}
	timer_stop(ts1);

	// calculate order parameter

	kuramoto_order_param(N,n,h,r,NULL);

	// wrap oscillator phases to [-pi,pi) [if that's is what you want]
	//
	// phase_wrap(m,h);

	// generate signal from phases

	for (size_t j=0; j<m; ++j) {
		x[j] = sin(h[j]);
	}

	// aggregate signal

	for (size_t k=0; k<n; ++k) {
		double yk = 0.0;
		for (size_t i=0; i<N; ++i) yk += x[N*k+i];
		y[k] = ooN*yk;
	}

	// write time stamp, order parameter and oscillator signals to file

	char ofile[] = "/tmp/kuramoto_audio.asc";     // output file (ASCII)
	printf("writing ASCII data to %s ...",ofile); fflush(stdout);
	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("Failed to open ASCII output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)(k+1)*dt); // time stamp
		fprintf(fp," %17.8f",r[k]);            // order parameter
		fprintf(fp," %17.8f",y[k]);            // aggregate signal
		for (size_t i=0; i<N; ++i) {
			fprintf(fp," %17.8f",x[N*k+i]);    // oscillator signal (waveform)
		}
		fprintf(fp,"\n");
	}
	if (fclose(fp) != 0) {
		perror("Failed to close ASCII output file");
		return EXIT_FAILURE;
	}
	printf(" done\n\n");

	// encode signals or aggregated signal as PCM and write to file
	//
	// if you have SoX, you can play the adio by, e.g.:
	//
	//   play -t raw -r 44.1k -e unsigned -b 16 -c 1 kuramoto_audio_44100Hz_c8a.u16
	//   play -t raw -r 44.1k -e float -b 32 -c 8 kuramoto_audio_44100Hz_c8.f32
	//
	// in the unaggregated case, you can play, e.g., just the 3rd channel with:
	//
	//   play -t raw -r 44.1k -e float -b 32 -c 8 kuramoto_audio_44100Hz_c8.f32 remix 3 0

	if (pcm) {
		const size_t smaxlen = 100;
		char rfile[smaxlen]; // raw PCM data file name
		if (pcm > 0) {
			snprintf(rfile,smaxlen,"/tmp/kuramoto_audio_%dHz_c%zu%s.u%d",F,N,cagg?"a":"",pcm);
			printf("writing PCM (%d-bit) data to %s ...",pcm,rfile); fflush(stdout);
		}
		else {
			snprintf(rfile,smaxlen,"/tmp/kuramoto_audio_%dHz_c%zu%s.f%d",F,N,cagg?"a":"",-pcm);
			printf("writing PCM (%d-bit floating-point) data to %s ...",-pcm,rfile); fflush(stdout);
		}
		FILE* const rp = fopen(rfile,"wb");
		if (rp == NULL) {
			perror("Failed to open PCM output file");
			return EXIT_FAILURE;
		}
		if (pcm_write(rp,cagg?y:x,cagg?n:N*n,pcm,-1.0,1.0) == -1) { // write PCM data
			return EXIT_FAILURE;
		}
		if (fclose(rp) != 0) {
			perror("Failed to close PCM output file");
			return EXIT_FAILURE;
		}
		printf(" done\n\n");
	}

	// if Gnuplot installed display order parameter and oscillator signals.
	// Else use your favourite plotting program on data in output file.

#ifdef _HAVE_GNUPLOT
	if (strncmp(gpterm,"noplot",7) != 0) {
		const size_t ns  = (size_t)ceil(Ts*f);  // start time step
		const size_t ne  = (size_t)ceil(Te*f);  // start time step
		char gfile[] = "/tmp/kuramoto_audio.gp"; // Gnuplot command file
		FILE* const gp = fopen(gfile,"w");
		if (gp == NULL) {
			perror("failed to open Gnuplot command file\n");
			return EXIT_FAILURE;
		}
		fprintf(gp,"set term \"%s\" title \"Kuramoto oscillator demo\" size 1600,1200\n",gpterm);
		fprintf(gp,"set xlabel \"time (secs)\"\n");
		fprintf(gp,"set ylabel \"mean phase\"\n");
		fprintf(gp,"set key right bottom Left rev\n");
		fprintf(gp,"# set grid\n");
		fprintf(gp,"set xr [%g:%g]\n",Ts,Te);
		fprintf(gp,"set yr [0:1.05]\n");
		fprintf(gp,"set ytics 0.5\n");
		fprintf(gp,"set multiplot layout 3,1\n");
		fprintf(gp,"set title \"Order parameter\"\n");
		fprintf(gp,"plot \"%s\" every ::%zu::%zu u 1:2 w l not\n",ofile,ns,ne);
		fprintf(gp,"set yr [-1.05:1.05]\n");
		fprintf(gp,"set title \"Aggregate oscillator signal (waveform)\"\n");
		fprintf(gp,"plot \"%s\" every ::%zu::%zu u 1:3 w l not\n",ofile,ns,ne);
		fprintf(gp,"set title \"Oscillator signals (waveforms)\"\n");
		fprintf(gp,"set ylabel \"amplitude\"\n");
		fprintf(gp,"plot \\\n");
		for (size_t i=0; i<N; ++i) fprintf(gp,"\"%s\" every ::%zu::%zu u 1:%zu w l not ,\\\n",ofile,ns,ne,4+i);
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

	free(y);
	free(x);
	free(r);
	free(h);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
