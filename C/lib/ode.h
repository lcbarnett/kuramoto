#ifndef ODE_H
#define ODE_H

// ODE solver macros; __VA_ARGS__ are the parameters to the ODE fun

#define ODE(ode,odefun,x,N,n,h,...) \
{ \
	switch (ode) { \
		case 1: { \
			printf("EULER : "#odefun); \
			double udot[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot,u,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i]; \
			}} \
			break; \
		case 2: { \
			printf("HEUN : "#odefun); \
			const double h2 = h/2.0; \
			double udot1[N], udot2[N]; \
			double v[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot1,u,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i]; \
				odefun(udot2,v,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]); \
			}} \
			break; \
		case 3: { \
			printf("RK4 : "#odefun); \
			const double h2 = h/2.0; \
			const double h6 = h/6.0; \
			double udot1[N],udot2[N],udot3[N],udot4[N]; \
			double v[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot1,u,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i]; \
				odefun(udot2,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i]; \
				odefun(udot3,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i]; \
				odefun(udot4,v,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]); \
			}} \
			break; \
	} \
}

// Chaotic systems

static inline void rossler(double* const xdot, const double* const x, const size_t N, const double a, const double b, const double c)
{
	xdot[0] = -(x[1]+x[2]);
	xdot[1] =  x[0]+a*x[1];
	xdot[2] =  b+x[2]*(x[0]-c);
}

static inline void lorenz(double* const xdot, const double* const x, const size_t N, const double s, const double r, const double b)
{
	xdot[0] = s*(x[1]-x[0]);
	xdot[1] = x[0]*(r-x[2])-x[1];
	xdot[2] = x[0]*x[1]-b*x[2];
}

static inline void thomas(double* const xdot, const double* const x, const size_t N, const double b)
{
	xdot[0] = -b*x[0] + sin(x[1]);
	xdot[1] = -b*x[1] + sin(x[2]);
	xdot[2] = -b*x[2] + sin(x[0]);
}

static inline void kmoto(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K)
{
	for (size_t i=0; i<N; ++i) {
		const double xi = x[i];
		const double* const Ki = K+N*i;
		double xdoti = w[i];
		for (size_t j=0; j<N; ++j) xdoti += Ki[j]*sin(x[j]-xi);
		xdot[i] = xdoti;
	}
}

static inline void kmotopl(double* const xdot, const double* const x, const size_t N, const double* const w, const double* const K, const double* const a)
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

static inline void lrnz96(double* const xdot, const double* const x, const size_t N, const double F)
{
	xdot[0] = (x[1]-x[N-2])*x[N-1]-x[0]+F;
	xdot[1] = (x[2]-x[N-1])*x[0]-x[1]+F;
	for (size_t i=2; i<N-1; ++i) xdot[i] = (x[i+1]-x[i-2])*x[i-1]-x[i]+F;
	xdot[N-1] = (x[0]-x[N-3])*x[N-2]-x[N-1]+F;
}

#endif // ODE_H
