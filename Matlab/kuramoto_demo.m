
% Default parameters (override on command line - see 'defvar.h')
%
% Note: phases are in radians on [0,2*pi), frequencies are in Hz

defvar('N',     4       ); % number of oscillators
defvar('Te',    100     ); % simulation equilibriation time - may be zero
defvar('T',     200     ); % simulation display time
defvar('dt',    0.01    ); % integration time increment
defvar('wmean', 0       ); % oscillator frequencies mean
defvar('wsdev', 0.2     ); % oscillator frequencies std. dev.
defvar('wseed', []      ); % oscillator frequencies random seed (empty for no seeding)
defvar('Kmean', 0.1     ); % oscillator coupling constants mean
defvar('Ksdev', Kmean/5 ); % oscillator coupling constants std. dev.
defvar('Kseed', []      ); % oscillator coupling constants random seed (empty for no seeding)
defvar('a',     []      ); % oscillator phase lag constant
defvar('hseed', []      ); % oscillator initial phases random seed (empty for no seeding)
defvar('nmean', 0.01    ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev', nmean/5 ); % oscillator input noise magnitude std. dev.
defvar('nseed', []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []      ); % oscillator input noise random seed (empty for no seeding)
defvar('unwrp', true    ); % oscillator phase plot: unwrapped and detrended? (else display on cylinder)
defvar('anim',  true    ); % order parameter animation?

% Times

ne = round(Te/dt);
n = round(T/dt);
assert(n > 0,'Simulation time too short, or time increment too large!');
T = n*dt; % adjusted simulation time

% Oscillator frequencies

if ~isempty(wseed), rstate = rng(wseed); end
w = wmean+wsdev*randn(N,1);
if ~isempty(wseed), rng(rstate); end
fprintf('\nOscillator natural frequencies =\n\n');
disp(w)

% Oscillator coupling constants, scaled by number of oscillators

if ~isempty(Kseed), rstate = rng(Kseed); end
K = (Kmean+Ksdev*randn(N))/N;
K(1:N+1:N^2) = 0;  % zero-out self-connections
if ~isempty(Kseed), rng(rstate); end
fprintf('\nOscillator coupling constants =\n\n');
disp(K)

% Uncorrelated Gaussian white noise, magnitudes log-normally distributed

if nmean > 0 % with input noise
	if ~isempty(nseed), rstate = rng(nseed); end
	lnv = log(1+nsdev^2/nmean^2);
	nmag = (lognrnd(log(nmean)-lnv/2,sqrt(lnv),1,N)); % normalise by sampling rate
	if ~isempty(nseed), rng(rstate); end
	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = nmag.*randn(ne+n,N); % convert to radians/sec
	if ~isempty(Iseed), rng(rstate); end
else
	I = []; % no input
end

% Initial phases, uniformly distributed on [0,2*pi)

if ~isempty(hseed), rstate = rng(hseed); end
h0 = 2*pi*rand(1,N);
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
r1   = r1  (ne+1:end)';
psi1 = psi1(ne+1:end)';

h2   = h2(:,ne+1:end)';
r2   = r2  (ne+1:end)';
psi2 = psi2(ne+1:end)';

t = linspace(0,T,n)';

% Display order parameters

figure(1); clf;
plot(t,[r1 r2]);
legend({'Euler','RK4'});
xlim([t(1) t(end)]);
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);
title(sprintf('\nKuramoto system: N = %d - order parameter magnitudes (r)\n',N));

% Display oscillator phases

figure(2); clf;
subplot(3,1,1);
if unwrp % Unwrapped, detrended
	h = ldetrend(h1,t);
	plot(t,h);
	title(sprintf('\nKuramoto system: N = %d - oscillator phases (unwrapped/detrended)\n',N));
	xlabel('time');
	ylabel('phase');
else      % On cylinder
	cylinder_plot(t,h1);
	title(sprintf('\nKuramoto system: N = %d - oscillator phases\n',N));
end
xlim([t(1) t(end)]);

% Display oscillator agregate signal

x = sin(h1);
y = mean(x,2);

subplot(3,1,2);
plot(t,y);
title(sprintf('\nKuramoto system: N = %d - oscillator agregate signal\n',N));
xlabel('time');
ylabel('magnitude');
xlim([t(1) t(end)]);
ylim([-1.05,+1.05]);

% Display all oscillator signals

subplot(3,1,3);
plot(t,x);
title(sprintf('\nKuramoto system: N = %d - oscillator signals\n',N));
xlabel('time');
ylabel('magnitude');
xlim([t(1) t(end)]);
ylim([-1.05,+1.05]);

% Display oscillator PSDs
%{
[psd,f] = pwelch(x);
subplot(3,1,3);
semilogy(f,psd);
title(sprintf('\nKuramoto system: N = %d - power spectral density\n',N));
xline(2*pi*dt*w);
xlim([0,pi]);
xticks([0,pi/4,pi/2,3*pi/4,pi]);
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
xlabel('angular frequency (rad/sec)');
ylabel('power (dB)');
%}

% Optionally display complex order parameter (animation)

if anim
	figure(3); clf;
	title(sprintf('\nKuramoto system: N = %d - complex order parameters (z)\n',N));
	clock_plot(t,[r1 r2],[psi1 psi2],{'Euler','RK4'},{'b','r'})
end
