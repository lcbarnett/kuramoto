#include <math.h>     // for maths functions
#include <string.h>   // for memcpy
#include <stdio.h>    // for perror

#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto.h" // for Kuramoto ODE solvers

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const size_t        n   = (size_t)mxGetScalar(prhs[0]); // number of integration increments
	const double        dt  =  mxGetScalar(prhs[1]);        // integration increment
	const double        b   =  mxGetScalar(prhs[2]);        // b parameter
	const double* const x0  =  mxIsEmpty(prhs[3]) ? NULL :  mxGetDoubles(prhs[3]);    // initial values
	const double* const I   =  mxIsEmpty(prhs[4]) ? NULL :  mxGetDoubles(prhs[4]);    // input (noise, etc.)
	const int           ode =  mxIsEmpty(prhs[5]) ? 1    : (int)mxGetScalar(prhs[5]); // 1 - Euler, 2 - Heun, 3 - RK4

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(3,n,mxREAL));

	// Initialise with input (noise, etc.) if present

	if (I) memcpy(x,I,3*n*sizeof(double));

	// Initial values if present

	if (x0) memcpy(x,x0,3*sizeof(double));

	// Run Rossler ODE solver

	switch (ode) {
		case 1: thomas_euler (n,dt,b,x); break;
		case 2: thomas_heun  (n,dt,b,x); break;
		case 3: thomas_rk4   (n,dt,b,x); break;
	}
}
