#ifndef KURAMATO_H
#define KURAMATO_H

static inline void kmoto_fun(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K)
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

static inline void kmotopl_fun(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K, const double* const a)
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

static inline void rossler_fun(double* const xdot, const double* const x, const double a, const double b, const double c)
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
	double* const  x  // the 3D variable (apprpriately initialised with noise or other input)
);

void rossler_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void rossler_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  x  // the 3D variable (apprpriately initialised with noise or other input)
);

static inline void lorenz_fun(double* const xdot, const double* const x, const double s, const double r, const double b)
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
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void lorenz_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void lorenz_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n,  // number of integration steps
	const   double h,  // integration increment
	const   double s,  // sigma parameter
	const   double r,  // rho   parameter
	const   double b,  // beta  parameter
	double* const  x   // the 3D variable (appropriately initialised with noise or other input)
);

static inline void thomas_fun(double* const xdot, const double* const x, const double b)
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
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void thomas_heun
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b  parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void thomas_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double b, // b parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

static inline void lrnz96_fun(double* const xdot, const double* const x, const size_t N, const double F)
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
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void lrnz96_heun
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

void lrnz96_rk4
(
	const   size_t N, // number of variables
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double F, // forcing parameter
	double* const  x  // the 3D variable (appropriately initialised with noise or other input)
);

// Utilities

static inline double phasewrap(const double x, const double u) // wrap to [-u,u)
{
	return x > 0.0 ? fmod(x+u,2.0*u)-u : fmod(x-u,2.0*u)+u;
}

void phase_wrap(const size_t m, double* const h, const double u); // wrap to [-u,u)

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
);

#endif // KURAMATO_H
