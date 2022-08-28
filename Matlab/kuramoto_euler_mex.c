#include "matrix.h"   // for Matlab matrix stuff
#include "kuramoto.h" // for Kuramoto ODE solvers

#define UNUSED __attribute__ ((unused))

// NOTE: C is row-major; since Matlab is column-major, transpose K matrix before calling the mex function

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	// read input parameters

	const size_t        N   = (size_t)*mxGetDoubles(prhs[0]);  // number of oscillators
	const size_t        n   = (size_t)*mxGetDoubles(prhs[1]);  // number of integration increments
	const double* const w   =  mxGetDoubles(prhs[2]);          // dt*frequencies
	const double* const K   =  mxGetDoubles(prhs[3]);          // dt*(coupling constants)
	const double        a   = *mxGetDoubles(prhs[4]);          // phase-lag (scalar)
	const double* const h0  =  mxGetDoubles(prhs[5]);          // initial oscillator phases

	// allocate output

	double* const h = mxGetDoubles(plhs[0] = mxCreateDoubleMatrix(N,n,mxREAL));

	// Euler method

	kuramoto_euler(N,n,w,K,a,h0,h);
}
