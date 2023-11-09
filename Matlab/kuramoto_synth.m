global seqno CD_SRATE MIDDLE_C

% Default parameters (override on command line - see 'defvar.h')

defvar('N',     8            ); % number of oscillators
defvar('T',     5            ); % sample time (seconds)
defvar('Td',    0.1          ); % sample display time (seconds)
defvar('fs',    CD_SRATE     ); % sampling frequency (Hz)
defvar('wmin',  0            ); % oscillator frequencies minimum (Hz)
defvar('wmax',  3*MIDDLE_C   ); % oscillator frequencies maximum (Hz)
defvar('wseed', []           ); % oscillator frequencies random seed (empty for no seeding)
defvar('kmean', 0.1          ); % oscillator coupling constants mean      (frequency multiplier - dimensionless)
defvar('ksdev', kmean/2      ); % oscillator coupling constants std. dev. (frequency multiplier - dimensionless)
defvar('kinhp', 0.1          ); % oscillator coupling probability of inhibitory connection
defvar('kseed', []           ); % oscillator couplinghttps://www.mathworks.com/help/releases/R2023a/matlab/ref/audiowrite.html?doclanguage=en&nocookie=true&prodfilter=ML%20SL%205G%20AE%20AT%20AA%20AU%20DR%20AS%20BI%20BL%20C2%20CM%20VP%20CT%20CF%20DB%20DF%20DD%20DH%20NN%20HS%20DS%20ET%20EC%20FH%20IT%20FI%20PO%20FL%20GD%20GC%20HD%20ES%20IA%20IP%20OT%20IC%20LP%20LS%20MG%20ME%20CO%20MJ%20MR%20TE%20DX%20AM%20MP%20MT%20NV%20OP%20DM%20PD%20AR%20PW%20PM%20RA%20RL%20RQ%20RB%20RP%20RF%20RK%20RO%20RC%20RR%20SI%20TF%20SX%20SQ%20SG%20SB%20SE%20SS%20BT%20LD%20PS%20SH%20MS%20VR%20VV%20CI%20RT%20SK%20SD%20CV%20SO%20DV%20PL%20XP%20SR%20SZ%20HW%20SF%20ST%20SM%20ZC%20ID%20TA%20UV%20VE%20VN%20VT%20WA%20LH%20WB%20WL&docviewer=helpbrowser&docrelease=R2023a&s_cid=pl_webdoc&loginurl=https%3A%2F%2F127.0.0.1%3A31515%2Ftoolbox%2Fmatlab%2Flogin%2Fweb%2Findex.html%3Fsnc%3DIANGSD%26external%3Dtrue%26channel%3D__mlfpmc__&searchsource=mw&snc=6TZA19&container=jshelpbrowser#d124e12665 constants random seed (empty for no seeding)
defvar('hseed', []           ); % oscillator initial phases random seed (empty for no seeding)
defvar('nmean', 0.01         ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev', nmean/5      ); % oscillator input noise magnitude std. dev.
defvar('nseed', []           ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []           ); % oscillator input noise random seed (empty for no seeding)
defvar('smode', 'Euler'      ); % simulation mode: 'Euler' or 'RK4'
defvar('wwin',  []           ); % Welch PSD window size (if empty set to number of samples/50)
defvar('codec', 'flac'       ); % audio codec
defvar('afseq', []           ); % audio file sequence number (empty to increment)
defvar('play',  false        ); % play audio?

assert(2*(N/2) == N,'Must be an even number of oscillators');

if isempty(afseq), seqno = seqno+1; else, seqno = afseq; end

% Times

n = round(T*fs);
assert(n > 0,'Simulation time too short!');
T = n/fs; % adjusted simulation time

nd = round(Td*fs);
assert(nd > 0,'Display time too short!');
Td = nd/fs; % adjusted display time

% Oscillator frequencies uniform random on [wmin,wmax)

if ~isempty(wseed), rstate = rng(wseed); end
w = wmin+(wmax-wmin)*rand(N,1);
if ~isempty(wseed), rng(rstate); end
fprintf('\nOscillator natural frequencies (Hz) =\n\n');
disp(w)

% Oscillator coupling multipliers are lognormal-distributed with given mean/standard deviation,
% inhibitory with given probability, scaled by number of oscillators

if ~isempty(kseed), rstate = rng(kseed); end
kln = log(1+ksdev^2/kmean^2);
K = w.*(lognrnd(log(kmean)-kln/2,sqrt(kln),N,N))/N;
K(1:N+1:N^2) = 0;  % zero-out self-connections
kinhib = rand(N)<kinhp;
K(kinhib) = -K(kinhib);
if ~isempty(kseed), rng(rstate); end
fprintf('Oscillator coupling constants =\n\n');
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

fprintf('Simulation: %s method ...',smode);
st = tic;
[h,r] = kuramoto(N,n,1/fs,w,K,[],h0,I,smode);
et = toc(st); fprintf(' %g seconds\n',et);

h = h';
r = r';

% Signal

x = sin(h);

left = 1:N/2;
right = N/2+1:N;
y = [mean(x(:,left), 2) mean(x(:,right),2)]; % left/right aggregate signal

% truncated for display

td = linspace(0,Td,nd)';
hd = h(1:nd,:);
rd = r(1:nd  );
xd = x(1:nd,:);
yd = y(1:nd,:);

% Display order parameters

figure(1); clf;
sgtitle(sprintf('Kuramoto synth: %d oscillators',N));

subplot(2,2,1);
plot(td,rd);
xlim([0 Td]);
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);
title(sprintf('\norder parameter magnitudes\n'));

% Display all oscillator signals

subplot(2,2,3);
plot(td,xd);
title(sprintf('\nOscillator signals\n'));
xlabel('time');
ylabel('magnitude');
xlim([0 Td]);
ylim([-1.05,+1.05]);

% Display aggregate signal

subplot(2,2,2);
plot(td,yd);
title(sprintf('\nStereo signal\n'));
xlabel('time');
ylabel('magnitude');
xlim([0 Td]);
ylim([-1.05,+1.05]);

% Display aggregate signal PSD

subplot(2,2,4);
if isempty(wwin), wwin = n/50; end
[psd,f] = pwelch(y,wwin,[],[],fs);
semilogy(f,psd);
title(sprintf('\nPower spectral density\n'));
xline(w);
xlim([0,2*wmax]);
xlabel('frequency (Hz)');
ylabel('power (dB)');

% Encode aggregate signal

if ~isempty(codec)
	afid = sprintf('kurasynth_mono_%03d.%s',seqno,codec);
	audiofile = fullfile(tempdir,afid);
	fprintf('\nWriting audio data to %s ...',audiofile);
	audiowrite(audiofile,y,fs);
	fprintf(' done\n\n');
	if play
		sound(y,fs);
	end
end

function [a,b] = betamv2ab(m,v)

	g = ((m*(1-m))/v)-1;
	a = g*m;
	b = g*(1-m);

end
