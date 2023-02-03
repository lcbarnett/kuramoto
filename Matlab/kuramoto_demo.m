
% Default parameters (override on command line - see 'defvar.h')

defvar('N',     20      ); % number of oscillators
defvar('wmean', 0       ); % oscillator frequencies mean
defvar('wsdev', pi/7    ); % oscillator frequencies std. dev.
defvar('wseed', []      ); % oscillator frequencies random seed (empty for no seeding)
defvar('Kmean', 0.8     ); % oscillator coupling constants mean
defvar('Ksdev', Kmean/6 ); % oscillator coupling constants std. dev.
defvar('Kseed', []      ); % oscillator coupling constants random seed (empty for no seeding)
defvar('a',     []      ); % oscillator phase lag constant
defvar('hseed', []      ); % oscillator initial phases random seed (empty for no seeding)
defvar('T',     200     ); % simulation time
defvar('dt',    0.01    ); % integration time increment
defvar('nmean', 0.05    ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev', nmean/3 ); % oscillator input noise magnitude std. dev.
defvar('nseed', []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []      ); % oscillator input noise random seed (empty for no seeding)
defvar('anim',  true    ); % order parameter animation?

% Random Kuramoto parameters

if ~isempty(wseed), rstate = rng(wseed); end
w = wmean + wsdev*randn(N,1); % oscillator frequencies normally distributed with mean wmean and std. dev wsdev
if ~isempty(wseed), rng(rstate); end

if ~isempty(Kseed), rstate = rng(Kseed); end
K = Kmean + Ksdev*randn(N); % coupling constants normally distributed with mean Kmean and std. dev. Ksdev
if ~isempty(Kseed), rng(rstate); end
K = K/N;                    % scale by system size

n = round(T/dt);
assert(n > 0,'Simulation time too short, or time increment too large!');
T = n*dt; % adjusted simulation time

if nmean > 0 % with input noise
	if ~isempty(nseed), rstate = rng(nseed); end
	lnv = log(1+nsdev^2/nmean^2);
	nmag = lognrnd(log(nmean)-lnv/2,sqrt(lnv),N,1); % per-oscillator noise magnitudes drawn from log-normal distribution
	if ~isempty(nseed), rng(rstate); end
end

% Generate input

if nmean > 0 % with input noise
	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = nmag.*randn(N,n); % uncorrelated Gaussian white noise
	if ~isempty(Iseed), rng(rstate); end
else
	I = zeros(N,n); % no input
end
if ~isempty(hseed), rstate = rng(hseed); end
I(:,1) = pi*(2*rand(N,1)-1);      % initial phases uniform on [-pi,pi]
if ~isempty(hseed), rng(rstate); end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('\n');

st1 = tic;
[h1,r1,psi1] = kuramoto(N,w,K,a,n,dt,I,'Euler');
et1 = toc(st1);
fprintf('Euler method : %g seconds\n',et1);

st2 = tic;
[h2,r2,psi2] = kuramoto(N,w,K,a,n,dt,I,'RK4');
et2 = toc(st2);
fprintf('RK4 method   : %g seconds\n',et2);

% Transpose stuff, etc., for plots

t   = linspace(0,T,n)';
r   = [r1;r2]';
psi = [psi1;psi2]';
h1  = h1';
h2  = h2';

% Display order parameter magnitude

figure(1); clf;
plot(t,r);
legend({'Euler','RK4'});
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);
title(sprintf('\nKuramoto system: N = %d - order parameter magnitudes (r)\n',N));

% Display oscillator phases on cylinder

figure(2); clf;
subplot(2,1,1);
cylinder_plot(t,h1);
title(sprintf('Euler\n'));
subplot(2,1,2);
cylinder_plot(t,h2);
title(sprintf('RK4\n'));
sgtitle(sprintf('\nKuramoto system: N = %d - oscillator phases\n',N));

% Optionally display complex order parameter (animation)

if anim
	figure(3); clf;
	title(sprintf('\nKuramoto system: N = %d - complex order parameters (z)\n',N));
	clock_plot(t,r,psi,{'Euler','RK4'},{'b','r'})
end
