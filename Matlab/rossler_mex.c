#include <math.h>     // for maths functions
#include <string.h>   // for memcpy
#include <stdio.h>    // for perror

#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto.h" // for Kuramoto ODE solvers

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const size_t        n  = (size_t)mxGetScalar(prhs[0]); // number of integration increments
	const double        dt =  mxGetScalar(prhs[1]);        // integration increment
	const double        a  =  mxGetScalar(prhs[2]);        // a parameter
	const double        b  =  mxGetScalar(prhs[3]);        // b parameter
	const double        c  =  mxGetScalar(prhs[4]);        // c parameter
	const double* const x0 =  mxIsEmpty(prhs[5]) ? NULL : mxGetDoubles(prhs[5]); // initial values
	const double* const I  =  mxIsEmpty(prhs[6]) ? NULL : mxGetDoubles(prhs[6]); // input (noise, etc.)

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(3,n,mxREAL));

	// Initialise with input (noise, etc.) if present

	if (I) memcpy(x,I,3*n*sizeof(double));

	// Initial values if present

	if (x0) memcpy(x,x0,3*sizeof(double));

	// Run Rossler ODE solver

	rossler_euler(n,dt,a,b,c,x);
}
