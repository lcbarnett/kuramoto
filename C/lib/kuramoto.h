#ifndef KURAMATO_H
#define KURAMATO_H

#define TWOPI (2.0*M_PI)

// Stuart-Landau: expected amplitude of sum of N oscillators
// uniform random on unit circle is approx SLMAGIC*sqrt(N)

#define SLMAGIC (0.887)

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

inline void kmoto_fun(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K)
{
	for (size_t i=0; i<N; ++i) {
		const double xi = x[i];
		const double* const Ki = K+N*i;
		double xdoti = w[i];
		for (size_t j=0; j<N; ++j) xdoti += Ki[j]*sin(x[j]-xi);
		xdot[i] = xdoti;
	}
}

void kmoto_euler
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	double* const         x  // oscillator phases, initialised with input
);

void kmoto_rk4
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	double* const         x  // oscillator phases, initialised with input
);

inline void kmotopl_fun(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K, const double* const a)
{
	for (size_t i=0; i<N; ++i) {
		const double xi = x[i];
		const double* const Ki = K+N*i;
		const double* const ai = a+N*i;
		double xdoti = w[i];
		for (size_t j=0; j<N; ++j) xdoti += Ki[j]*sin(x[j]-xi-ai[j]);
		xdot[i] = xdoti;
	}
}

void kmotopl_euler
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	const   double* const a, // phase lags
	double* const         x  // oscillator phases, initialised with input
);

void kmotopl_rk4
(
	const   size_t        N, // number of oscillators
	const   size_t        n, // number of integration increments
	const   double        h, // integration increment
	const   double* const w, // frequencies
	const   double* const K, // coupling constants
	const   double* const a, // phase lags
	double* const         x  // oscillator phases, initialised with input
);

inline void rossler_fun(double* const xdot, const double* const x, const double a, const double b, const double c)
{
	xdot[0] = -(x[1]+x[2]);
	xdot[1] =  x[0]+a*x[1];
	xdot[2] =  b+x[2]*(x[0]-c);
}

void rossler_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  u  // the 3D variable (apprpriately initialised with noise or other input)
);

void rossler_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  u  // the 3D variable (apprpriately initialised with noise or other input)
);

inline void lorenz_fun(double* const xdot, const double* const x, const double s, const double r, const double b)
{
	xdot[0] = s*(x[1]-x[0]);
	xdot[1] = x[0]*(r-x[2])-x[1];
	xdot[2] = x[0]*x[1]-b*x[2];
}

void lorenz_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
);

void lorenz_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n,  // number of integration steps
	const   double h,  // integration increment
	const   double s,  // sigma parameter
	const   double r,  // rho   parameter
	const   double b,  // beta  parameter
	double* const  u   // the 3D variable (appropriately initialised with noise or other input)
);

inline void thomas_fun(double* const xdot, const double* const x, const double b)
{
	xdot[0] = -b*x[0] + sin(x[1]);
	xdot[1] = -b*x[1] + sin(x[2]);
	xdot[2] = -b*x[2] + sin(x[0]);
}

void thomas_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
);

void thomas_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
);

inline void lrnz96_fun(double* const xdot, const double* const x, const size_t N, const double F)
{
	xdot[0] = (x[1]-x[N-2])*x[N-1]-x[0]+F;
	xdot[1] = (x[2]-x[N-1])*x[0]-x[1]+F;
	for (size_t i=2; i<N-1; ++i) xdot[i] = (x[i+1]-x[i-2])*x[i-1]-x[i]+F;
	xdot[N-1] = (x[0]-x[N-3])*x[N-2]-x[N-1]+F;
}

void lrnz96_euler
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
);

void lrnz96_rk4
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
);

// Utilities

static inline double phasewrap(const double x, const double u) // wrap to [-u,u)
{
	return x > 0.0 ? fmod(x+u,2.0*u)-u : fmod(x-u,2.0*u)+u;
}

void phase_wrap(const size_t m, double* const h, const double u); // wrap to [-u,u)

#endif // KURAMATO_H
