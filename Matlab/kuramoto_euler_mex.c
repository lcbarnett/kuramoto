#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto.h" // for Kuramoto ODE solvers
#include <string.h>   // for memcpy

#define UNUSED __attribute__ ((unused))

// NOTE: C is row-major; since Matlab is column-major, transpose K matrix before calling the mex function

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const size_t        N   = (size_t)mxGetScalar(prhs[0]);  // number of oscillators
	const size_t        n   = (size_t)mxGetScalar(prhs[1]);  // number of integration increments
	const double* const w   =  mxGetDoubles(prhs[2]);        // dt*frequencies
	const double* const K   =  mxGetDoubles(prhs[3]);        // dt*(coupling constants)
	const double        a   =  mxGetScalar(prhs[4]);         // phase-lag (scalar)
	const double* const h0  =  mxGetDoubles(prhs[5]);        // initial oscillator phases
	const double* const I   =  mxIsEmpty(prhs[6]) ? NULL : mxGetDoubles(prhs[6]); // sqrt(dt)*input

	// allocate output

	double* const h = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(N,n,mxREAL));

	// Initialise with input (if present)
	//
	// Note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes (this is what we want!)

	if (I) { // have input
		memcpy(h,I,N*n*sizeof(double)); // initialise phases with input
	}

	// Initial phases at t = 0 (plus input!)

	for (size_t i=0; i<N; ++i) {
		h[i] += h0[i];
	}

	// Euler method

	kuramoto_euler(N,n,w,K,a,h);
}
