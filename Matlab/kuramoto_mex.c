#include <math.h>     // for maths functions
#include <string.h>   // for memcpy

#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto.h" // for Kuramoto ODE solvers

#define UNUSED __attribute__ ((unused))

// NOTE: C is row-major; since Matlab is column-major, transpose K matrix before calling the mex function

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const   size_t N   = (size_t)mxGetScalar(prhs[0]); // number of oscillators
	const   size_t n   = (size_t)mxGetScalar(prhs[1]); // number of integration increments
	const   double dt  =  mxGetScalar(prhs[2]);        // time increment (secs)
	double* const  w   =  mxGetDoubles(prhs[3]);       // frequencies (radians/sec)
	double* const  K   =  mxGetDoubles(prhs[4]);       // coupling constants (radians/sec)
	double* const  a   =  mxIsEmpty(prhs[5]) ? NULL : mxGetDoubles(prhs[5]);     // phase-lags (radians)
	double* const  I   =  mxIsEmpty(prhs[6]) ? NULL : mxGetDoubles(prhs[6]);     // input (radians)
	const   int    RK4 =  mxIsEmpty(prhs[7]) ? 1    : (int)mxGetScalar(prhs[7]); // use RK4? (default: yes)

	// allocate output

	double* const h = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(N,n,mxREAL));

	// Initialise with input (if present)
	//
	// Note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes (this is what we want!)

	if (I) { // have input
		memcpy(h,I,N*n*sizeof(double)); // initialise phases with input
	}

	// run Kuramoto ODE solver

	if (RK4) {
		double* const k1 = mxCalloc(4*N,sizeof(double));
		if (a) kuramoto_rk4pl(N,n,dt,w,K,a,h,k1); else kuramoto_rk4(N,n,dt,w,K,h,k1);
		mxFree(k1);
	}
	else {
		if (a) kuramoto_eulerpl(N,n,dt,w,K,a,h); else kuramoto_euler(N,n,dt,w,K,h);
	}
}
