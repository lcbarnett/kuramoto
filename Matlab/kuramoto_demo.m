
% Default parameters (override on command line - see 'defvar.h')

defvar('N',     10      ); % number of oscillators
defvar('Te',    20      ); % simulation equilibriation time - may be zero (secs)
defvar('T',     200     ); % simulation display time (secs)
defvar('dt',    0.01    ); % integration time increment (secs)
defvar('fs',    1000    ); % sampling frequency (arbitrary, Hz)
defvar('wmean', 200     ); % oscillator frequencies mean (Hz)
defvar('wsdev', wmean/2 ); % oscillator frequencies std. dev. (Hz)
defvar('wseed', []      ); % oscillator frequencies random seed (empty for no seeding)
defvar('Kmean', 100     ); % oscillator coupling constants mean (Hz)
defvar('Ksdev', Kmean/2 ); % oscillator coupling constants std. dev. (Hz)
defvar('Kiprob',0.1     ); % oscillator coupling constants inhibitory connection probability
defvar('Kseed', []      ); % oscillator coupling constants random seed (empty for no seeding)
defvar('a',     []      ); % oscillator phase lag constant
defvar('hseed', []      ); % oscillator initial phases random seed (empty for no seeding)
defvar('nmean', 10      ); % oscillator input noise magnitude mean (Hz, zero for no noise)
defvar('nsdev', nmean/4 ); % oscillator input noise magnitude std. dev. (Hz)
defvar('nseed', []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []      ); % oscillator input noise random seed (empty for no seeding)
defvar('unwrp', true    ); % oscillator phase plot: unwrapped and detrended? (else display on cylinder)
defvar('anim',  true    ); % order parameter animation?

% Times

ne = round(Te/dt);
n = round(T/dt);
assert(n > 0,'Simulation time too short, or time increment too large!');
T = n*dt; % adjusted simulation time

% Oscillator frequencies, beta-distributed on (0,fs)

wm = wmean/fs; % normalise to [0,1)
ws = wsdev/fs; % normalise to [0,1)
wg = (wm*(1-wm))/(ws^2)-1;
assert(wg > 0,'Frequencies variance too big!');
if ~isempty(wseed), rstate = rng(wseed); end
w = 2*pi*betarnd(wm*wg,(1-wm)*wg,N,1); % convert to radians/sec
if ~isempty(wseed), rng(rstate); end
fprintf('\nOscillator natural frequencies =\n');
disp(fs*w/2/pi)

% Oscillator coupling constants, beta-distributed on (0,fs), some inhibitory, scaled by number of oscillators

km = Kmean/fs; % normalise to [0,1)
ks = Ksdev/fs; % normalise to [0,1)
kg = (km*(1-km))/(ks^2)-1;
if ~isempty(Kseed), rstate = rng(Kseed); end
K = 2*pi*betarnd(km*kg,(1-km)*kg,N)/N; % convert to radians/sec
i = find(rand(N)<Kiprob);              % connections inhibitory with specified probability
K(i) = -K(i);                          % reverse signs of inhibitory connections
K(1:N+1:N^2) = 0;                      % zero-out self-connections
if ~isempty(Kseed), rng(rstate); end
fprintf('Inhibitory connections: %d of %d\n',length(i),N*(N-1));

% Uncorrelated Gaussian white noise, magnitudes log-normally distributed

if nmean > 0 % with input noise
	if ~isempty(nseed), rstate = rng(nseed); end
	lnv = log(1+nsdev^2/nmean^2);
	nmag = (lognrnd(log(nmean)-lnv/2,sqrt(lnv),1,N))/fs; % normalise by sampling rate
	if ~isempty(nseed), rng(rstate); end

	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = 2*pi*nmag.*randn(ne+n,N); % convert to radians/sec
	if ~isempty(Iseed), rng(rstate); end
else
	I = []; % no input
end

 % Initial phases, uniformly distributed on [-pi,pi]

if ~isempty(hseed), rstate = rng(hseed); end
h0 = pi*(2*rand(1,N)-1);
if ~isempty(hseed), rng(rstate); end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('\n');

st1 = tic;
[h1,r1,psi1] = kuramoto(N,ne+n,dt,w,K,a,h0,I,'Euler');
et1 = toc(st1);
fprintf('Euler method : %g seconds\n',et1);

st2 = tic;
[h2,r2,psi2] = kuramoto(N,ne+n,dt,w,K,a,h0,I,'RK4');
et2 = toc(st2);
fprintf('RK4 method   : %g seconds\n',et2);

% Truncate equilibriation and transpose stuff, etc., for plots

h1   = h1(:,ne+1:end)';
h2   = h2(:,ne+1:end)';
r1   = r1(ne+1:end)';
r2   = r2(ne+1:end)';
psi1 = psi1(ne+1:end)';
psi2 = psi2(ne+1:end)';

t = linspace(0,T,n)';

% Display order parameters

figure(1); clf;
plot(t,[r1 r2]);
legend({'Euler','RK4'});
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);
title(sprintf('\nKuramoto system: N = %d - order parameter magnitudes (r)\n',N));

% Display oscillator phases

figure(2); clf;
if unwrp

	% Detrended

	x1 = detrend(h1);
	x2 = detrend(h2);
	subplot(2,1,1);
	plot(t,x1);
	title(sprintf('Euler\n'));
	subplot(2,1,2);
	plot(t,x2);
	title(sprintf('RK4\n'));
	sgtitle(sprintf('\nKuramoto system: N = %d - oscillator phases (unwrapped/detrended)\n',N));

else

	% On cylinder

	figure(2); clf;
	subplot(2,1,1);
	cylinder_plot(t,h1);
	title(sprintf('Euler\n'));
	subplot(2,1,2);
	cylinder_plot(t,h2);
	title(sprintf('RK4\n'));
	sgtitle(sprintf('\nKuramoto system: N = %d - oscillator phases\n',N));

end

% Optionally display complex order parameter (animation)

if anim
	figure(3); clf;
	title(sprintf('\nKuramoto system: N = %d - complex order parameters (z)\n',N));
	clock_plot(t,[r1 r2],[psi1 psi2],{'Euler','RK4'},{'b','r'})
end
