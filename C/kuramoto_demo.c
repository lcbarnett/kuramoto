#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "kuramoto.h"

// random double in range

inline double randd(const double a, const double b)
{
	return a + (b-a)*((double)rand()/(double)RAND_MAX);
}

// Main function

int main()
// int main(int argc, char* argv[])
{
	const size_t N  = 4;
	const double T  = 20.0;
	const double dt = 0.01;
	const size_t n  = (size_t)round(T/dt);

	double* const w = calloc(N,  sizeof(double));
	double* const K = calloc(N*N,sizeof(double));
	double* const h = calloc(N*n,sizeof(double));

	// set up some frequencies

	const double wrange = M_PI/7.0;
	for (size_t i=0; i<N; ++i) {
		w[i] = randd(-wrange,wrange);
	}

	// set up some connectivity

	const double Kmean = 0.8/(double)N;
	const double Krange = Kmean/6.0;
	for (size_t i=0; i<N; ++i) {
		for (size_t j=0; j<N; ++j) {
			if (i != j) {
				K[N*i+j] = randd(Kmean-Krange,Kmean+Krange);
			}
		}
	}

	// set up some input noise

	const double Irange = M_PI/20.0;
	for (size_t k=0; k<N*n; ++k) {
		h[k] = randd(-Irange,Irange);
	}

	kuramoto_euler(N,n,w,K,h);

	// write data to file

	FILE* const fp = fopen("/tmp/kuramto_demo.out","w");
	if (fp == NULL) {
		perror("Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		fprintf(fp,"%17.8f",(double)k*dt);
		for (size_t i=0; i<N; ++i) {
			fprintf(fp,"\t%17.8f",h[n*i+k]);
		}
		fprintf(fp,"\n");
	}

	if (fclose(fp) != 0) {
		perror("Failed to close output file");
		return EXIT_FAILURE;
	}

	free(h);
	free(K);
	free(w);

	return EXIT_SUCCESS;
}
