#include <math.h>   // for maths functions
#include <stdlib.h> // for malloc, etc.
#include <stdio.h>  // for malloc, etc.

#include "kuramoto.h"

// NOTE 1: C is row-major; bear in mind when writing interfaces! E.g. for
// Matlab (column-major) you should transpose the matrices K and a before calling.
//
// NOTE 2: All phases (including h) are dimensionless; multiply by 2π for radians.
//
// NOTE 3: To initialise input phase variables with Wiener noise, scale by sqrt(dt).
// For other (deterministic) input scale by dt.

void kuramoto_euler	// Euler method
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double hti = ht[i];
			double dhti = wdt[i];
			for (size_t j=0; j<N; ++j) dhti += Kdti[j]*sin(ht[j]-hti);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_eulerpl // Euler method with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	const   double* const a,   // phase lags (dimensionless)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double dhti = wdt[i];
			for (size_t j=0; j<N; ++j) dhti += Kdti[j]*sin(ht[j]-hti-ai[j]);
			ht1[i] += hti+dhti; // update next time step (adding in input already in ht1)
		}
	}
}

void kuramoto_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t  N,          // number of oscillators
	const   size_t  n,          // number of integration increments
	const   double* const  wdt, // frequencies x dt (dimensionless)
	const   double* const  Kdt, // coupling constants x dt (dimensionless)
	double* const   k1,         // buffer for RK4 coefficients (size must be 4*N)
	double* const   h           // oscillator phases, initialised with input (dimensionless)
)
{
	// set up coefficients buffers

	double* const k2 = k1+N;
	double* const k3 = k2+N;
	double* const k4 = k3+N;

	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double hti = ht[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]-hti);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k1[j]-htk1i);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k2[j]-htk2i);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k3[j]-htk3i);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
	}
}

void kuramoto_rk4pl // Classic Runge-Kutta (RK4) with phase lags
(
	const   size_t        N,   // number of oscillators
	const   size_t        n,   // number of integration increments
	const   double* const wdt, // frequencies x dt (dimensionless)
	const   double* const Kdt, // coupling constants x dt (dimensionless)
	const   double* const a,   // phase lags (dimensionless)
	double* const         k1,  // buffer for RK4 coefficients (size must be 4*N)
	double* const         h    // oscillator phases, initialised with input (dimensionless)
)
{
	double* const k2 = k1+N;
	double* const k3 = k2+N;
	double* const k4 = k3+N;

	// ODE solver

	for (double* ht=h; ht<h+N*(n-1); ht+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double hti = ht[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]-hti-ai[j]);
			k1[i] = ki;
		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk1i = ht[i]+k1[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k1[j]-htk1i-ai[j]);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk2i = ht[i]+k2[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k2[j]-htk2i-ai[j]);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Kdti = Kdt+N*i;
			const double* const ai = a+N*i;
			const double htk3i = ht[i]+k3[i];
			double ki = wdt[i];
			for (size_t j=0; j<N; ++j) ki += Kdti[j]*sin(ht[j]+k3[j]-htk3i-ai[j]);
			k4[i] = ki;
		}

		// update next time step (adding in input already in ht1)
		double* const ht1 = ht+N;
		for (size_t i=0; i<N; ++i) ht1[i] += ht[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
	}
}

void kuramoto_order_param // calculate order parameter magnitude/phase
(
	const   size_t         N, // number of oscillators
	const   size_t         n, // number of integration increments
	const   double* const  h, // oscillator phases
	double* const          r, // order parameter magnitude
	double* const          p  // order parameter phase (NULL if not required)
)
{
	const double OON = 1.0/(double)N;
	double* rt = r;
	double* pt = p;
	for (const double* ht=h; ht<h+N*n; ht+=N) {
		double xmt = 0.0;
		double ymt = 0.0;
		for (size_t i=0; i<N; ++i) {
#ifdef _GNU_SOURCE
			double c,s;
			sincos(ht[i],&s,&c);
			xmt += c;
			ymt += s;
#else
			xmt += cos(ht[i]);
			ymt += sin(ht[i]);
#endif
		}
		xmt *= OON;
		ymt *= OON;
		*rt++ = hypot(xmt,ymt);
		if (p) *pt++ = atan2(ymt,xmt);
	}
}

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
)
{
	// scale parameters appropriately according to time increment

	const double sqrtdt = sqrt(dt);
	for (size_t i=0; i<N;   ++i) w[i] *= dt;
	for (size_t i=0; i<N*N; ++i) K[i] *= dt;
	for (size_t i=0; i<N;   ++i) a[i] *= dt;
	for (size_t i=N; i<N*n; ++i) x[i] *= sqrtdt; // cf. Ornstein-Uhlenbeck process
	for (size_t i=N; i<N*n; ++i) y[i] *= sqrtdt; // cf. Ornstein-Uhlenbeck process

	// adjust growth constants by mean coupling

	for (size_t i=0; i<N; ++i) {
		const double* const Ki = K+N*i;
		double Kbari = 0.0;
		for (size_t j=0; j<N; ++j) Kbari += Ki[j];
		a[i] -= Kbari;
	}

	// ODE solver

	for (double *xt=x, *yt = y; xt<x+N*(n-1); xt+=N, yt+=N) {
		double* const xt1 = xt+N;
		double* const yt1 = yt+N;
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xti*xti+yti*yti);
			double dxti = vi*xti - w[i]*yti;
			double dyti = vi*yti + w[i]*xti;
			for (size_t j=0; j<N; ++j) dxti += Ki[j]*xt[j];
			for (size_t j=0; j<N; ++j) dyti += Ki[j]*yt[j];
			xt1[i] += xti+dxti; // update next time step (adding in input already in yt1)
			yt1[i] += yti+dyti; // update next time step (adding in input already in xt1)
		}
	}
}

void stulan_magnitudes // calculate magnitudes of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r  // oscillator magnitude
)
{
	double* rt = r;
	for (double *xt=x, *yt = y; xt<x+N*n; ++xt, ++yt) *rt++ = hypot(*xt,*yt);
}

void stulan_phases // calculate phases of oscillators
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  h  // oscillator phase
)
{
	double* ht = h;
	for (double *xt=x, *yt = y; xt<x+N*n; ++xt, ++yt) *ht++ = atan2(*yt,*xt);
}

void stulan_order_param // calculate order parameter magnitude/phase
(
	const   size_t N, // number of oscillators
	const   size_t n, // number of integration increments
	double* const  x, // oscillator real part
	double* const  y, // oscillator imag part
	double* const  r, // order parameter magnitude
	double* const  p  // order parameter phase (NULL if not required)
)
{
	const double OON = 1.0/(double)N;
	double* rt = r;
	double* pt = p;
	for (double *xt=x, *yt = y; xt<x+N*n; xt+=N, yt+=N) {
		double xmt = 0.0;
		for (size_t i=0; i<N; ++i)  xmt += xt[i];
		xmt *= OON;
		double ymt = 0.0;
		for (size_t i=0; i<N; ++i)  ymt += yt[i];
		ymt *= OON;
		*rt++ = hypot(xmt,ymt);
		if (p) *pt++ = atan2(ymt,xmt);
	}
}

/*
void stulan_rk4 // Classic Runge-Kutta (RK4)
(
	const size_t        N,  // number of oscillators
	const size_t        n,  // number of integration increments
	const double        dt, // time integration step
	const double* const w,  // dt*frequencies
	const double* const K,  // dt*frequencies*(coupling constants)/N
	const double* const a,  // dt*(growth constants) - (K summed over 2nd index)
	double*       const x,  // oscillator real part, to be computed by numerical ODE (pre-initialised with input)
	double*       const y,  // oscillator imag part, to be computed by numerical ODE (pre-initialised with input)
	double*       const k   // buffer for RK4 coefficients (size must be 8*N)
)
{
	double* const kx1 = k  +N;
	double* const ky1 = kx1+N;
	double* const kx2 = ky1+N;
	double* const ky2 = kx2+N;
	double* const kx3 = ky2+N;
	double* const ky3 = kx3+N;
	double* const kx4 = ky3+N;
	double* const ky4 = kx4+N;

	double* yt=y;
	for (double* xt=x; xt<x+N*(n-1); xt+=N,yt+=N) {

		// k1
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xt[i]*xt[i]+yt[i]*yt[i]);
			double kxi = vi*xti - w[i]*yt[i];
			double kyi = vi*yti + w[i]*xt[i];
			for (size_t j=0; j<N; ++j) kxi += Ki[j]*xt[j];
			for (size_t j=0; j<N; ++j) kyi += Ki[j]*yt[j];
			kx1[i] = kxi;
			ky1[i] = kyi;

		}

		// k2
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double xti = xt[i];
			const double yti = yt[i];
			const double vi = a[i]-dt*(xt[i]*xt[i]+yt[i]*yt[i]);

			const double xtk1i = xt[i]+kx1[i];
			double kxi = vi*xti - w[i]*yt[i];
			double kyi = vi*yti + w[i]*xt[i];
			for (size_t j=0; j<N; ++j) kxi += Ki[j]*(xt[j]+kx1[j]-xtk1i);
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k1[j]-htk1i);
			k2[i] = ki/2.0;
		}

		// k3
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk2i = xt[i]+k2[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k2[j]-htk2i);
			k3[i] = ki/2.0;
		}

		// k4
		for (size_t i=0; i<N; ++i) {
			const double* const Ki = K+N*i;
			const double htk3i = xt[i]+k3[i];
			double ki = w[i];
			for (size_t j=0; j<N; ++j) ki += Ki[j]*sin(xt[j]+k3[j]-htk3i);
			k4[i] = ki;
		}

		// update next time step (adding in input already in xt1)
		double* const xt1 = xt+N;
		for (size_t i=0; i<N; ++i) {
			xt1[i] += xt[i] + (k1[i]+4.0*k2[i]+4.0*k3[i]+k4[i])/6.0;
		}
	}
}
*/

void rossler_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  u  // the 3D variable (apprpriately initialised with noise or other input)
)
{
	for (double* v=u; v<u+3*(n-1); v+=3) {
		v[3] += v[0] - h*(v[1]+v[2]);
		v[4] += v[1] + h*(v[0]+a*v[1]);
		v[5] += v[2] + h*(b+v[2]*(v[0]-c));
	}
}

void rossler_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double a, // a parameter
	const   double b, // b parameter
	const   double c, // c parameter
	double* const  u  // the 3D variable (apprpriately initialised with noise or other input)
)
{
	const double hh = h/2.0;
	const double h6 = h/6.0;

	double k10,k11,k12,k20,k21,k22,k30,k31,k32,k40,k41,k42;

	for (double* v=u; v<u+3*(n-1); v+=3) {

		k10 = -(v[1]+v[2]);
		k11 =  v[0]+a*v[1];
		k12 =  b+v[2]*(v[0]-c);

		k20 = -((v[1]+hh*k11)+(v[2]+hh*k12));
		k21 =  (v[0]+hh*k10)+a*(v[1]+hh*k11);
		k22 =  b+(v[2]+hh*k12)*((v[0]+hh*k10)-c);

		k30 = -((v[1]+hh*k21)+(v[2]+hh*k22));
		k31 =  (v[0]+hh*k20)+a*(v[1]+hh*k21);
		k32 =  b+(v[2]+hh*k22)*((v[0]+hh*k20)-c);

		k40 = -((v[1]+h*k31)+(v[2]+h*k32));
		k41 =  (v[0]+h*k30)+a*(v[1]+h*k31);
		k42 =  b+(v[2]+h*k32)*((v[0]+h*k30)-c);

		v[3] += v[0] + h6*(k10+2.0*k20+2.0*k30+k40);
		v[4] += v[1] + h6*(k11+2.0*k21+2.0*k31+k41);
		v[5] += v[2] + h6*(k12+2.0*k22+2.0*k32+k42);
	}
}

void lorenz_euler
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
)
{
	for (double* v=u; v<u+3*(n-1); v+=3) {
		v[3] += v[0] + h*(s*(v[1]-v[0]));
		v[4] += v[1] + h*(v[0]*(r-v[2])-v[1]);
		v[5] += v[2] + h*(v[0]*v[1]-b*v[2]);
	}
}

void lorenz_rk4 // Classic Runge-Kutta (RK4)
(
	const   size_t n, // number of integration steps
	const   double h, // integration increment
	const   double s, // sigma parameter
	const   double r, // rho   parameter
	const   double b, // beta  parameter
	double* const  u  // the 3D variable (appropriately initialised with noise or other input)
)
{
	const double hh = h/2.0;
	const double h6 = h/6.0;

	double k10,k11,k12,k20,k21,k22,k30,k31,k32,k40,k41,k42;

	for (double* v=u; v<u+3*(n-1); v+=3) {

		k10 = s*(v[1]-v[0]);
		k11 = v[0]*(r-v[2])-v[1];
		k12 = v[0]*v[1]-b*v[2];

		k20 = s*((v[1]+hh*k11)-(v[0]+hh*k10));
		k21 = (v[0]+hh*k10)*(r-(v[2]+hh*k12))-(v[1]+hh*k11);
		k22 = (v[0]+hh*k10)*(v[1]+hh*k11)-b*(v[2]+hh*k12);

		k30 = s*((v[1]+hh*k21)-(v[0]+hh*k20));
		k31 = (v[0]+hh*k20)*(r-(v[2]+hh*k22))-(v[1]+hh*k21);
		k32 = (v[0]+hh*k20)*(v[1]+hh*k21)-b*(v[2]+hh*k22);

		k40 = s*((v[1]+h*k31)-(v[0]+h*k30));
		k41 = (v[0]+h*k30)*(r-(v[2]+h*k32))-(v[1]+h*k31);
		k42 = (v[0]+h*k30)*(v[1]+h*k31)-b*(v[2]+h*k32);

		v[3] += v[0] + h6*(k10+2.0*k20+2.0*k30+k40);
		v[4] += v[1] + h6*(k11+2.0*k21+2.0*k31+k41);
		v[5] += v[2] + h6*(k12+2.0*k22+2.0*k32+k42);
	}
}

// Utilities

void phase_wrap(const size_t m, double* const h, const double u) // wrap to [-u,u)
{
	for (size_t k=0; k<m; ++k) h[k] = phasewrap(h[k],u);
}
