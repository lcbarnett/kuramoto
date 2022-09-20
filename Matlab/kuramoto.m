function [h,r,psi,T,n] = kuramoto(N,w,K,a,h0,T,dt,V)

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
% T     simulation time                       (positive double; or, if negative, number of integration time steps is n = -T)
% dt    integration time increment            (positive double)
% V     mode/input noise variance/covariance  ('Euler', 'RK4', or positive scalar, vector or square positive-definite matrix of size N)
%
% h     oscillator phases (unwrapped)         (N x n matrix)
% r     order parameter magnitude             (row vector of length n)
% psi   order parameter phase (wrapped)       (row vector of length n)
% T     simulation time (possibly adjusted)   (positive double)
% n     integration time steps                (positive integer)
%
% NOTE 1: K(i,j) is connection strength from oscillator j to oscillator i.
%
% NOTE 2: If T is entered as a negative integer, then -T is taken as the number of time increments.
%
% NOTE 3: Euler method is faster (by a factor of about 5), but RK4 is more accurate.
%
% NOTE 4: Noisy method is Euler with Gaussian noise input (no point using RK4 with noise, so not implemented!)
%
% NOTE 5: To wrap the oscillator phases h to [-pi,pi), do:
%
%     h = mod(h+pi,2*pi)-pi;
%
% See kuramoto_demo.m script for usage.
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

assert(isa(dt,'double') && isscalar(dt) && dt > 0,'Integration increment must be a positive scalar double');

T = double(T);
assert(isscalar(T),'Integration simulation time must be scalar');
if T > eps
	n = floor(T/dt);
	assert(n > 0,'Simulation time too short, or time increment too large!');
elseif T < -eps
	n = -T;
	assert(floor(n) == n && n > 0,'Number of integration steps must be a positive scalar integer');
else
	error('Bad simulation time');
end
T = n*dt; % adjusted simulation time (<= T)

% Call mex ODE simulation (returned phase matrix h is N x n)
%
% We transpose K so that K(i,j) is connection strength j --> i

if ischar(V)
	switch upper(V)
		case 'RK4',   h = kuramoto_rk4_mex(N,n,w*dt,K'*dt,a,h0);
		case 'EULER', h = kuramoto_euler_mex(N,n,w*dt,K'*dt,a,h0);
		otherwise,    error('Unknown simulation mode');
	end
else
	assert(isa(V,'double'),'Noise variance/covariance must be a positive scalar double, a vector of positive doubles, or a positive-definite matrix of doubles matching the specified number of oscillators');
	if isscalar(V)
		assert(V >= 0,'Noise variance must be positive');
		L = sqrt(V)*ones(N); % for uncorrelated white noise
	elseif isvector(V)
		assert(all(V >= 0),'Noise variances must be positive');
		L = diag(sqrt(V));
	else
		assert(ismatrix(V) && size(V,1) == N && size(V,2) == N,'Noise variance/covariance must be a positive scalar double, a vector of positive doubles, or a positive-definite matrix of doubles matching the specified number of oscillators');
		[L,cholp] = chol(V,'lower');
		assert(cholp == 0,'Noise variance/covariance matrix must be positive-definite');
	end
	I = L*randn(N,n); % Input: Gaussian white noise
	h = kuramoto_noisy_mex(N,n,w*dt,K'*dt,a,h0,I*sqrt(dt));
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
