function x = lorenz(n,dt,parms,x0,I,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for Lorenz C mex function - all input checking happens here!
%
% To compile "mex" files, see Makefile in this directory
%
% n        number of time increments
% dt       time integration step
% parms    Lorenz sigma, rho, beta parameters (3-vector)
% x0       initial values (3-vector, or empty)
% I        input noise (n x 3 matrix, or empty)
% mode     simulation mode: 'Euler', 'Heun', or 'RK4'
%
% x        Rossler 3D variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input checks - essential, as the mex functions do NO error checking!

n = double(n);
assert(isscalar(n) && floor(n) == n && n > 0,'Number of time increments must be a positive scalar integer');

assert(isa(dt,'double') && isscalar(dt) && dt > 0,'Integration increment must be a positive scalar double');

if nargin < 3 || isempty(parms)
	parms = [10, 28, 8/3];
else
	assert(isa(parms,'double') && isvector(parms) && length(parms) == 3,'Rossler parameters must be empty (for defaults), or a 3-vector of doubles');
end

if nargin < 4 || isempty(x0)
	x0 = [1,1,1];
else
	assert(isempty(x0) || (isa(x0,'double') && isvector(x0) && length(x0) == 3),'Initial values must be empty, or a 3-vector of doubles');
end

if nargin < 5
	I = [];
else
	assert(isempty(I) || (isa(I,'double') && ismatrix(I) && size(I,1) == n && size(I,2) == 3),'Input noise must be empty, or an n x 3 matrix of doubles');
end

if nargin < 6 || isempty(mode)
	ode = 2;
else
	assert(ischar(mode),'Simulation mode must be empty, ''Euler'', ''Heun'', or ''RK4''');
	switch upper(mode)
		case 'EULER', ode = 1;
		case 'HEUN',  ode = 2;
		case 'RK4',   ode = 3;
		otherwise,    error('Simulation mode must be empty, ''Euler'', ''Heun'', or ''RK4''');
	end
end

% Call mex ODE simulation (returned matrix x is n x 3)
%
% Input noise is Wiener (Brownian), so scaled by sqrt(dt); cf. Ito simulation of Ornstein-Uhlenbeck process

x = lorenz_mex(n,dt,parms(1),parms(2),parms(3),x0,I*sqrt(dt),ode);
