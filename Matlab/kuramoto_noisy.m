function [h,r,psi,T,n] = kuramoto_noisy(N,w,K,a,V,h0,T,dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for noisy kuramoto C mex function - all input checking happens here!
%
% To compile "mex" files, see Makefile in this directory
%
% N     number of oscillators                (positive integer)
% w     oscillator frequencies               (vector of length N)
% K     oscillator coupling constants        (scalar or square matrix of size N)
% a     phase lag                            (scalar)
% V     input noise variance/covariance      (positive scalar or square positive-definite matrix of size N)
% h0    initial phases of oscillators        (scalar or vector of length N)
% T     simulation time                      (positive double; or, if negative, number of integration time steps is n = -T)
% dt    integration time increment           (positive double)
%
% h     oscillator phases (unwrapped)        (N x n matrix)
% r     order parameter magnitude            (row vector of length n)
% psi   order parameter phase (wrapped)      (row vector of length n)
% T     simulation time (possibly adjusted)  (positive double)
% n     integration time steps               (positive integer)
%
% NOTE 1: K(i,j) is connection strength from oscillator j to oscillator i.
%
% NOTE 2: Euler simulation with noise (Ito stochastic ODE-style); since there is
%         noise, we don't really need the extra accuracy of RK4, etc. Here we use
%         Gaussian noise, but could be adapted for other noise processes.
%
% NOTE 3: To wrap the oscillator phases h to [-pi,pi), do:
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

assert(isa(V,'double'),'Variance/covariance must be a positive scalar double, or a positive-definite matrix of doubles matching the specified number of oscillators');
if isscalar(V)
	assert(V > 0,'Variance must be positive');
	L = sqrt(V)*ones(N); % for uncorrelated white noise
else
	assert(ismatrix(V) && size(V,1) == N && size(V,2) == N,'Variance/covariance must be a positive scalar double, or a positive-definite matrix of doubles matching the specified number of oscillators');
	[L,cholp] = chol(V,'lower');
	assert(cholp == 0,'Variance/covariance matrix must be positive-definite');
end

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

I = L*randn(N,n); % Input: Gaussian white noise

% Call mex stochastic ODE Ito simulation (returned phase matrix h is N x n);
% note that noise is scaled by sqrt(dt) to gett the variance scaling right
%
% We transpose K so that K(i,j) is connection strength j --> i

h = kuramoto_noisy_mex(N,n,w*dt,K'*dt,a,I*sqrt(dt),h0);

% Order parameter (if requested)

if nargout > 1
	x = mean(cos(h));
	y = mean(sin(h));
	r = hypot(x,y);
	if nargout > 2
		psi = atan2(y,x); % wrapped on [-pi,pi)
	end
end
