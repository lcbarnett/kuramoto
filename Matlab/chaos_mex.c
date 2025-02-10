#include <math.h>     // for maths functions
#include <string.h>   // for memcpy
#include <stdio.h>    // for perror

#include "matrix.h"   // for Matlab matrix stuff
#include "ode.h"      // for chaotic ODE solvers

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const int           sys = (int)mxGetScalar(prhs[1]);    // 1 - Lorenz, 2 - Rossler, 3 - Thomas
	const int           ode = (int)mxGetScalar(prhs[2]);    // 1 - Euler,  2 - Heun,    3 - RK4
	const size_t        n   = (size_t)mxGetScalar(prhs[3]); // number of integration increments
	const double        dt  =  mxGetScalar(prhs[4]);        // integration increment
	const double* const p   =  mxGetDoubles(prhs[5]);       // parameters
	const double* const x0  =  mxIsEmpty(prhs[6]) ? NULL :  mxGetDoubles(prhs[6]);  // initial values
	const double* const I   =  mxIsEmpty(prhs[7]) ? NULL :  mxGetDoubles(prhs[7]);  // input (noise, etc.)

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(3,n,mxREAL));

	// Initialise with input (noise, etc.) if present

	if (I) memcpy(x,I,3*n*sizeof(double));

	// Initial values if present

	if (x0) memcpy(x,x0,3*sizeof(double));

	// Run ODE solver

	switch (sys) {
		case 1: ODE(ode,lorenz, x,3,n,dt,p[0],p[1],p[2]); break;
		case 2: ODE(ode,rossler,x,3,n,dt,p[0],p[1],p[2]); break;
		case 3: ODE(ode,thomas, x,3,n,dt,p[0]);           break;
	}

}
