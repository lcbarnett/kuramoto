#include <math.h>     // for maths functions
#include <string.h>   // for memcpy

#include "matrix.h"   // for Matlab matrix stuff
#include "mex.h"   // for Matlab matrix stuff

#include "ode.h"      // for chaotic ODE solvers

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const sys_t         sys = (sys_t)mxGetScalar(prhs[0]);  // LORENZ = 0, ROSSLER = 1, THOMAS = 2
	const ode_t         ode = (ode_t)mxGetScalar(prhs[1]);  // EULER  = 0, HEUN    = 1, RKFOUR = 2
	const size_t        n   = (size_t)mxGetScalar(prhs[2]); // number of integration increments
	const double        dt  =  mxGetScalar(prhs[3]);        // integration increment
	const double* const p   =  mxGetDoubles(prhs[4]);       // parameters
	const double* const x0  =  mxIsEmpty(prhs[5]) ? NULL :  mxGetDoubles(prhs[5]);  // initial values
	const double* const I   =  mxIsEmpty(prhs[6]) ? NULL :  mxGetDoubles(prhs[6]);  // input (noise, etc.)

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const x = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(3,n,mxREAL));

	// Initialise with input (noise, etc.) if present

	if (I) memcpy(x,I,3*n*sizeof(double));

	// Initial values if present

	if (x0) memcpy(x,x0,3*sizeof(double));

	// Run ODE solver
	switch (sys) {
		case LORENZ  : {ODE(ode,lorenz, x,3,n,dt,p[0],p[1],p[2]);} break;
		case ROSSLER : {ODE(ode,rossler,x,3,n,dt,p[0],p[1],p[2]);} break;
		case THOMAS  : {ODE(ode,thomas, x,3,n,dt,p[0]);}           break;
		default      : mexErrMsgIdAndTxt("chaos_mex:badsys","Unknown system type");
	}
}
