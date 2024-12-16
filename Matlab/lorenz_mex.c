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
	const double        s   =  mxGetScalar(prhs[2]);        // sigma parameter
	const double        r   =  mxGetScalar(prhs[3]);        // rho   parameter
	const double        b   =  mxGetScalar(prhs[4]);        // beta  parameter
	const double* const x0  =  mxIsEmpty(prhs[5]) ? NULL :  mxGetDoubles(prhs[5]);    // initial values
	const double* const I   =  mxIsEmpty(prhs[6]) ? NULL :  mxGetDoubles(prhs[6]);    // input (noise, etc.)
	const int           RK4 =  mxIsEmpty(prhs[7]) ? 1    : (int)mxGetScalar(prhs[7]); // use RK4? (else Euler)

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(3,n,mxREAL));

	// Initialise with input (noise, etc.) if present

	if (I) memcpy(x,I,3*n*sizeof(double));

	// Initial values if present

	if (x0) memcpy(x,x0,3*sizeof(double));

	// Run Lorenz ODE solver

	if (RK4) lorenz_rk4(n,dt,s,r,b,x); else lorenz_euler(n,dt,s,r,b,x);
}
