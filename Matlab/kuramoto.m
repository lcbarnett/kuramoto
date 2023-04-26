function [h,r,psi] = kuramoto(N,n,dt,w,K,a,h0,I,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for kuramoto C mex functions - all input checking happens here!
%
% To compile "mex" files, see Makefile in this directory
%
% N     number of oscillators                 (positive integer)
% n     number of time increments             (positive integer)
% dt    time integration step                 (scalar)                             : seconds
% w     oscillator frequencies                (scalar or vector of length N)       : radians/second
% K     oscillator coupling constants         (scalar or square matrix of size N)  : radians/second
% a     phase lags                            (scalar or square matrix of size N)  : radians
% h0    initial phases                        (vector of length N)                 : radians
% I     input noise                           (n x N matrix or empty for no input) : radians/sqrt(second)
% mode  simulation mode                       ('Euler' or 'RK4')                   : string
%
% h     oscillator phases (unwrapped)         (N x n matrix)                       : radians
% r     order parameter magnitude             (row vector of length n)             : dimensionless
% psi   order parameter phase (wrapped)       (row vector of length n)             : radians
%
% NOTE 1: Phases are interpreted as radians on [0,2*pi), and frequencies angular; i.e., in radians/second
%
% NOTE 2: K(i,j) is connection strength from oscillator j to oscillator i.
%
% NOTE 3: Euler method is faster (by a factor of about 5), but RK4 is more accurate.
%
% NOTE 4: To wrap the oscillator phases h to [-pi,pi), do:
%
%     h = mod(h+pi,2*pi)-pi;
%
% See kuramoto_demo.m script for example usage.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input checks - essential, as the mex functions do NO error checking!

N = double(N);
assert(isscalar(N) && floor(N) == N && N > 0,'Number of oscillators must be a positive scalar integer');

n = double(n);
assert(isscalar(n) && floor(n) == n && n > 0,'Number of time increments must be a positive scalar integer');

assert(isa(dt,'double') && isscalar(dt) && dt > 0,'Integration increment must be a positive scalar double');

assert(isa(w,'double'),'Frequencies must be a scalar double or a vector of doubles matching the specified number of oscillators');
if isscalar(w)
	w = w*ones(N,1);
else
	assert(isvector(w) && length(w) == N,'Frequencies must be a scalar double or a vector of doubles matching the specified number of oscillators');
end

assert(isa(K,'double'),'Coupling constants must be a scalar double, or a square matrix of doubles matching the specified number of oscillators');
if isscalar(K)
	K = K*ones(N);    % uniform coupling
	K(1:N+1:N*N) = 0; % zeros on diagonal!
else
	assert(ismatrix(K) && size(K,1) == N && size(K,2) == N,'Coupling constants must be a scalar double, or a square matrix of doubles matching the specified number of oscillators');
end

if ~isempty(a)
	assert(isa(a,'double'),'Phase lags must be empty, a scalar double, or a square matrix of doubles matching matching the specified number of oscillators');
	if isscalar(a)
		a = a*ones(N);    % uniform phase lag
		a(1:N+1:N*N) = 0; % zeros on diagonal!
	else
		assert(ismatrix(a) && size(a,1) == N && size(a,2) == N,'Phase lags must be empty, a scalar double, or a square matrix of doubles matching the specified number of oscillators');
	end
end

assert(isempty(h0) || (isa(h0,'double') && isvector(h0) && length(h0) == N),'Initial phases must be empty, or a vector of doubles of length N');

assert(isempty(I) || (isa(I,'double') && ismatrix(I) && size(I,1) == n && size(I,2) == N),'Input noise must be empty, or an n x N matrix of doubles');

if isempty(mode)
	RK4 = 1;
else
	assert(ischar(mode),'Simulation mode must be empty, or ''Euler'' or ''RK4''');
	switch upper(mode)
		case 'RK4',   RK4 = 1;
		case 'EULER', RK4 = 0;
		otherwise,    error('Unknown simulation mode; must be empty, ''Euler'' or ''RK4''');
	end
end

% Call mex ODE simulation (returned phase matrix h is n x N)
%
% We transpose K so that K(i,j) is connection strength j --> i, same for a
%
% Input noise is Weiner (Brownian), so scaled by sqrt(dt); cf. Ito simulation of Ornstein-Uhlenbeck process

h = kuramoto_mex(N,n,w*dt,K'*dt,a',h0,I*sqrt(dt),RK4);

% Order parameter (if requested)

if nargout > 1
	x = mean(cos(h));
	y = mean(sin(h));
	r = hypot(x,y);
	if nargout > 2
		psi = atan2(y,x); % wrapped on [-pi,pi)
	end
end
