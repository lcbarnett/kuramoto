#ifndef KURAMATO_H
#define KURAMATO_H

#define TWOPI (2.0*M_PI)

// Stuart-Landau: expected amplitude of sum of N oscillators
// uniform random on unit circle is approx SLMAGIC*sqrt(N)

#define SLMAGIC (0.887)

// Experimental /////////////////////////////////////////////////////////////////

typedef double* darray;

static inline darray* matalloc(size_t rows, size_t cols, const darray buffer) // allocate a (row-major indexed) matrix of doubles
{
	// If a buffer is supplied, it is *essential* that it is of length at least rows*cols
	darray* x = malloc(rows*sizeof(double*));
	if (x == NULL) {
		perror("memory allocation failed\n");
		return NULL;
	}
	if (buffer == NULL) { // allocate buffer - x must be deallocated by matfree(x)
		x[0] = calloc(rows*cols,sizeof(double)); // zero-initialises
		if (x[0] == NULL) {
			perror("memory allocation failed\n");
			return NULL;
		}
	}
	else {
		x[0] = buffer; // attach to supplied buffer - x must be deallocated by free(x)
	}
	for (size_t i=1; i<rows; ++i) x[i] = x[i-1] + cols;
	return x; // so x[i][j] is entry in i-th row, j-th column
}

static inline void matfree(darray* x)
{
	free(x[0]);
	free(x);
}

void kuramoto_euler_alt	// Euler method
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   darray* const Kdt, // coupling constants x dt (radians)
	darray* const         h    // oscillator phases, initialised with input (radians)
);

/////////////////////////////////////////////////////////////////////////////////

// Kuramoto model

void kuramoto_euler	// Euler method
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_eulerpl // Euler method with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	const   double* const a,   // phase lags (radians)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	double* const         k1,  // buffer for RK4 coefficients (size must be 4*N)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (radians)
	const   double* const Kdt, // coupling constants x dt (radians)
	const   double* const a,   // phase lags (radians)
	double* const         k1,  // buffer for RK4 coefficients (size must be 4*N)
	double* const         h    // oscillator phases, initialised with input (radians)
);

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
);

// Stuart-Landau model

void stulan_euler // Euler method
(
	const   size_t N,  // number of oscillators
	const   size_t n,  // number of integration increments
	const   double dt, // time increment (secs)
	double* const  w,  // frequencies (1/sec)
	double* const  K,  // coupling constants (1/sec)
	double* const  a,  // growth constants (1/sec)
	double* const  x,  // oscillator real part (dimensionless), initialised with input
	double* const  y   // oscillator imag part (dimensionless), initialised with input
);

void stulan_magnitudes // calculate magnitudes of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r  // oscillator magnitude
);

void stulan_phases // calculate phases of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  h  // oscillator phase
);

void stulan_order_param // calculate order parameter magnitude/phase
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r, // order parameter magnitude
	double* const  p  // order parameter phase (NULL if not required)
);

// Utilities

static inline double wmpi2pi(const double x) // wrap to [-pi,pi)
{
	return x > 0.0 ? fmod(x+M_PI,2.0*M_PI)-M_PI : fmod(x-M_PI,2.0*M_PI)+M_PI;
}

void phase_wrap(const size_t m, double* const h);

#endif // KURAMATO_H
