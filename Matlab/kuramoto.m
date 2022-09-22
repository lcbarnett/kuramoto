function [h,r,psi] = kuramoto(N,w,K,a,h0,n,dt,I,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for kuramoto C mex functions - all input checking happens here!
%
% To compile "mex" files, see Makefile in this directory
%
% N     number of oscillators                 (positive integer)
% w     oscillator frequencies                (vector of length N)
% K     oscillator coupling constants         (scalar or square matrix of size N)
% a     phase lag                             (scalar)
% h0    initial phases of oscillators         (scalar or vector of length N)
% n     number of time increments             (positive integer)
% dt    integration time increment            (positive double)
% I     input                                 (N x n matrix or empty for no input)
% mode  simulation mode                       ('Euler' or 'RK4')
%
% h     oscillator phases (unwrapped)         (N x n matrix)
% r     order parameter magnitude             (row vector of length n)
% psi   order parameter phase (wrapped)       (row vector of length n)
%
% NOTE 1: K(i,j) is connection strength from oscillator j to oscillator i.
%
% NOTE 2: Euler method is faster (by a factor of about 5), but RK4 is more accurate.
%
% NOTE 3: Input is scaled by sqrt(dt), as per Ito SDE simulation (cf. Ornstein-Uhlenbeck)
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

assert(isa(w,'double') && isvector(w) && length(w) == N,'Frequencies must be a vector of doubles matching the specified number of oscillators');

assert(isa(K,'double'),'Coupling constants must be a scalar double, or a vector of doubles matching the specified number of oscillators');
if isscalar(K)
	K = K*ones(N); % uniform coupling
else
	assert(ismatrix(K) && size(K,1) == N && size(K,2) == N,'Coupling constants must be a scalar double, or a square matrix of doubles matching matching the specified number of oscillators');
end

assert(isa(a,'double') && isscalar(a),'Phase lag must be scalar double');

assert(isa(h0,'double'),'Initial oscillator phases must be a vector of doubles matching the specified number of oscillators');
if isscalar(h0)
	h0 = h0*ones(N,1);
else
	assert(isvector(h0)  && length(h0) == N,'Initial oscillator phases must be a vector of doubles matching the specified number of oscillators');
end

n = double(n);
assert(isscalar(n) && floor(n) == n && n > 0,'Number of time increments must be a positive scalar integer');

assert(isa(dt,'double') && isscalar(dt) && dt > 0,'Integration increment must be a positive scalar double');

assert(isempty(I) || (isa(I,'double') && ismatrix(I) && size(I,1) == N && size(I,2) == n),'Input must be empty, or an N x n matrix of doubles');

% Call mex ODE simulation (returned phase matrix h is N x n)
%
% Note: we transpose K so that K(i,j) is connection strength j --> i

assert(ischar(mode),'Simulation mode must be ''Euler'' or ''RK4''');
switch upper(mode)
case 'RK4'
	h = kuramoto_rk4_mex(N,n,w*dt,K'*dt,a,h0,I*sqrt(dt));
case 'EULER'
	h = kuramoto_euler_mex(N,n,w*dt,K'*dt,a,h0,I*sqrt(dt));
otherwise
	error('Unknown simulation mode; must be ''Euler'' or ''RK4''');
end

% Order parameter (if requested)

if nargout > 1
	x = mean(cos(h));
	y = mean(sin(h));
	r = hypot(x,y);
	if nargout > 2
		psi = atan2(y,x); % wrapped on [-pi,pi)
	end
end
