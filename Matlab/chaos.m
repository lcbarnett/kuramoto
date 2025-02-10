function x = chaos(sys,ode,n,dt,parms,x0,I)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for Lorenz C mex function - all input checking happens here!
%
% To compile "mex" files, see Makefile in this directory
%
% sys     'Lorenz', 'Rossler', or 'Thomas'
% ode     'Euler', 'Heun', or 'RK4'
% n        number of time increments
% dt       time integration step
% parms    system parameters (3-vector)
% x0       initial values (3-vector, or empty)
% I        input noise (n x 3 matrix, or empty)
%
% x        Rossler 3D variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input checks - essential, as the mex functions do NO error checking!

assert(ischar(sys),'Simulation mode must be ''Lorenz'', ''Rossler'', or ''Thomas''');
switch upper(sys)
	case 'LORENZ',  sys = 1;
	case 'ROSSLER', sys = 2;
	case 'THOMAS',  sys = 3;
	otherwise,      error('Simulation mode must be ''Lorenz'', ''Rossler'', or ''Thomas''');
end

assert(ischar(ode),'Simulation mode must be ''Euler'', ''Heun'', or ''RK4''');
switch upper(ode)
	case 'EULER', ode = 1;
	case 'HEUN',  ode = 2;
	case 'RK4',   ode = 3;
	otherwise,    error('Simulation mode must be ''Euler'', ''Heun'', or ''RK4''');
end

n = double(n);
assert(isscalar(n) && floor(n) == n && n > 0,'Number of time increments must be a positive scalar integer');

assert(isa(dt,'double') && isscalar(dt) && dt > 0,'Integration increment must be a positive scalar double');

if nargin < 5 || isempty(parms)
	switch (sys)
		case 1, parms(1) = 10.00; parms(2) = 28.00; parms(3) = 8.0/3.0;
		case 2, parms(1) =  0.10; parms(2) =  0.10; parms(3) =   14.00;
		case 3, parms(1) =  0.20; parms(2) =   NaN; parms(3) =     NaN;
	end
else
	switch (sys)
		case 1, assert(isa(parms,'double') && isvector(parms) && length(parms) == 3,'Lorenz parameters must be empty (for defaults), or a 3-vector of doubles');
		case 2, assert(isa(parms,'double') && isvector(parms) && length(parms) == 3,'Rossler parameters must be empty (for defaults), or a 3-vector of doubles');
		case 3, assert(isa(parms,'double') && isvector(parms),'Thomas parameter must be empty (for default), or a scalar double'); parms(2) = NaN; parms(3) = NaN; % ignore last two parameters!
	end
end

if nargin < 6 || isempty(x0)
	x0 = [1,1,1];
else
	assert(isempty(x0) || (isa(x0,'double') && isvector(x0) && length(x0) == 3),'Initial values must be empty, or a 3-vector of doubles');
end

if nargin < 5
	I = [];
else
	assert(isempty(I) || (isa(I,'double') && ismatrix(I) && size(I,1) == n && size(I,2) == 3),'Input noise must be empty, or an n x 3 matrix of doubles');
end

% Call mex ODE simulation (returned matrix x is n x 3)
%
% Input noise is Wiener (Brownian), so scaled by sqrt(dt); cf. Ito simulation of Ornstein-Uhlenbeck process

x = chaos_mex(sys,ode,n,dt,parms,x0,I*sqrt(dt),ode);
