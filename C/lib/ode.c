#include <stddef.h>

#include "ode.h"

const char* ode2str(const ode_t ode)
{
	static char* odes[] = {"Euler", "Heun", "RK4"};
	if (ode < 0 || ode > RKFOUR) return NULL;
	return odes[ode];
}

const char* sys2str(const sys_t sys)
{
	if (sys < 0 || sys > LNINESIX) return NULL;
	static char* syss[] = {"Lorenz", "Rossler", "Thomas", "Kuramoto", "KuramotoPL", "Lorenz96"};
	return syss[sys];
}
