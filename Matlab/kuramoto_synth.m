CD_SRATE = 44100.0
MIDDLE_C = 261.625565

% Default parameters (override on command line - see 'defvar.h')

defvar('N',     8          ); % number of oscillators
defvar('T',     10         ); % sample time (seconds)
defvar('fs',    CD_SRATE   ); % sampling frequency (Hz)
defvar('wmin',  0          ); % oscillator frequencies minimum (Hz)
defvar('wmax',  3*MIDDLE_C ); % oscillator frequencies maximum (Hz)
defvar('wseed', []         ); % oscillator frequencies random seed (empty for no seeding)
defvar('kmean', 0.1        ); % oscillator coupling constants mean      (frequency multiplier - dimensionless)
defvar('ksdev', kmean/2    ); % oscillator coupling constants std. dev. (frequency multiplier - dimensionless)
defvar('kinhp', 0.1        ); % oscillator coupling probability of inhibitory connection
defvar('kseed', []         ); % oscillator coupling constants random seed (empty for no seeding)
defvar('hseed', []         ); % oscillator initial phases random seed (empty for no seeding)
defvar('nmean', 0.01       ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev', nmean/5    ); % oscillator input noise magnitude std. dev.
defvar('nseed', []         ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []         ); % oscillator input noise random seed (empty for no seeding)

% Times

n = round(T*fs);
assert(n > 0,'Simulation time too short, or time increment too large!');
T = n/fs; % adjusted simulation time

% Oscillator frequencies uniform random on [wmin,wmax)

if ~isempty(wseed), rstate = rng(wseed); end
w = wmin+(wmax-wmin)*rand(N,1);
if ~isempty(wseed), rng(rstate); end
fprintf('\nOscillator natural frequencies (Hz) =\n\n');
disp(w)

% Oscillator coupling multipliers are gamma-distributed with given mean/standard deviation,
% inhibitory with given probability, then scaled by number of oscillators

if ~isempty(kseed), rstate = rng(kseed); end
kln = log(1+ksdev^2/kmean^2);
K = w.*(lognrnd(log(kmean)-kln/2,sqrt(kln),N,N)); % normalise by sampling rate
K(1:N+1:N^2) = 0;  % zero-out self-connections
kinhib = rand(N)<kinhp;
K(kinhib) = -K(kinhib);
if ~isempty(kseed), rng(rstate); end
fprintf('\nOscillator coupling constants =\n\n');
disp(K)

% Uncorrelated Gaussian white noise, magnitudes log-normally distributed

if nmean > 0 % with input noise
	if ~isempty(nseed), rstate = rng(nseed); end
	nln = log(1+nsdev^2/nmean^2);
	nmag = (lognrnd(log(nmean)-nln/2,sqrt(nln),1,N));
	if ~isempty(nseed), rng(rstate); end
	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = nmag.*randn(n,N);
	if ~isempty(Iseed), rng(rstate); end
else
	I = []; % no input
end

% Initial phases, uniformly distributed on [0,2*pi)

if ~isempty(hseed), rstate = rng(hseed); end
h0 = 2*pi*rand(1,N);
if ~isempty(hseed), rng(rstate); end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('\n'); st = tic;
[h,r] = kuramoto(N,n,1/fs,w,K,[],h0,I,'Euler');
et = toc(st); fprintf('Euler method : %g seconds\n',et);

h   = h';
r   = r';

% Signal

x = sin(h);

t = linspace(0,T,n)';

% Display order parameters

figure(1); clf;
plot(t,r);
xlim([t(1) t(end)]);
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);
title(sprintf('\nKuramoto system: N = %d - order parameter magnitudes (r)\n',N));

return;

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

function [a,b] = betamv2ab(m,v)

	g = ((m*(1-m))/v)-1;
	a = g*m;
	b = g*(1-m);

end
