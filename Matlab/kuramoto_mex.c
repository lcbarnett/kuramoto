#include <math.h>     // for maths functions
#include <string.h>   // for memcpy
#include <stdio.h>    // for perror

#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto_old.h" // for Kuramoto ODE solvers

#define UNUSED __attribute__ ((unused))

// NOTE: C is row-major; since Matlab is column-major, transpose Kdt matrix before calling the mex function

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const   size_t        N    = (size_t)mxGetScalar(prhs[0]); // number of oscillators
	const   size_t        n    = (size_t)mxGetScalar(prhs[1]); // number of integration increments
	const   double* const wdt  =  mxGetDoubles(prhs[2]);       // frequencies x dt (radians)
	const   double* const Kdt  =  mxGetDoubles(prhs[3]);       // coupling constants x dt (radians)
	const   double* const a    =  mxIsEmpty(prhs[4]) ? NULL : mxGetDoubles(prhs[4]);     // phase-lags (radians)
	const   double* const h0   =  mxIsEmpty(prhs[5]) ? NULL : mxGetDoubles(prhs[5]);     // initial phases (radians)
	const   double* const Isdt =  mxIsEmpty(prhs[6]) ? NULL : mxGetDoubles(prhs[6]);     // process noise x sqrt(dt)  (radians)
	const   int           RK4  =  mxIsEmpty(prhs[7]) ? 1    : (int)mxGetScalar(prhs[7]); // use RK4? (default: yes)

	// allocate output (note that if compiling with -R2018a, mxCreateDoubleMatrix zero-initializes)

	double* const h = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(N,n,mxREAL));

	// Initialise phases with process noise (if present)

	if (Isdt) memcpy(h,Isdt,N*n*sizeof(double)); // initialise phases with noise

	// Initialise phases with initial values (if present)

	if (h0) memcpy(h,h0,N*sizeof(double));

	// run Kuramoto ODE solver

	if (RK4) {
		double* const k1 = mxCalloc(4*N,sizeof(double));
		if (a) kuramoto_rk4pl(N,n,wdt,Kdt,a,k1,h); else kuramoto_rk4(N,n,wdt,Kdt,k1,h);
		mxFree(k1);
	}
	else {
		if (a) kuramoto_eulerpl(N,n,wdt,Kdt,a,h); else kuramoto_euler(N,n,wdt,Kdt,h);
	}
}
