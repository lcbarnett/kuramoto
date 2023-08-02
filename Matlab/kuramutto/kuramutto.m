global seqno CD_SRATE

% Sampling rate

fs = CD_SRATE;

% Default parameters (override on command line - see 'defvar.m')

defvar('T',     5            ); % sample time (seconds)
defvar('Td',    0.1          ); % sample display time (seconds)
defvar('fspec', 'kmfreq.m'   ); % frequency specs file
defvar('cspec', 'kmconn.m'   ); % connectivity specs file
defvar('nspec', 'kmnois.m'   ); % nois specs file
defvar('hseed', []           ); % oscillator initial phases random seed (empty for no seeding)
defvar('smode', 'Euler'      ); % simulation mode: 'Euler' or 'RK4'
defvar('pdto',  1            ); % polynomial detrend order
defvar('codec', 'flac'       ); % audio codec
defvar('afseq', []           ); % audio file sequence number (empty to increment)
defvar('play',  true         ); % play audio?

if isempty(afseq), seqno = seqno+1; else, seqno = afseq; end

% Times

ndt = round(T*fs);
assert(ndt > 0,'Simulation time too short!');
T = ndt/fs; % adjusted simulation time
t = linspace(0,T,ndt)';

nddt = round(Td*fs);
assert(nddt > 0,'Display time too short!');
Td = nddt/fs; % adjusted display time
td = linspace(0,Td,nddt)';

% Construct the mutt

kmbody;

% Oscillator frequencies

kmfreq;

% Oscillator coupling constants

kmconn;

% Input noise

kmnois;

I = inoise;

% Initial phases, uniformly distributed on [0,2*pi)

if ~isempty(hseed), rstate = rng(hseed); end
h0 = 2*pi*rand(1,nosc);
if ~isempty(hseed), rng(rstate); end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('Simulation: %s method ...',smode);
st = tic;
[h,r] = kuramoto(nosc,ndt,1/fs,w,K,[],h0,I,smode);
et = toc(st); fprintf(' %g seconds\n',et);

h = h';
r = r';

% Voice

lhead = 1:3;
rhead = 4:6;
head = [lhead rhead];
nhead = length(head);
arse  = nosc-nhead+1:nosc;

fprintf('\nPolynomial detrend at order %d\n',pdto);
p = detrend(h(:,arse),pdto);
p = (p-mean(p))./std(p);
p = p./max(abs(p));

x = p.*sin(h(:,head)); % modulate with detrended phases :-)
x = x./max(abs(x));

y = [mean(x(:,lhead), 2) mean(x(:,rhead),2)]; % left/right aggregate signal
y = y./max(abs(y));

vfreq = 4;
venv = abs(sin(vfreq*2*pi*fs*t));
y = venv.*y;

% truncated for display

rd = r(1:nddt  );
pd = p(1:nddt,:);
xd = x(1:nddt,:);
yd = y(1:nddt,:);

% Display order parameters

figure(1); clf;
sgtitle('Kuro-mutt-o');

subplot(2,2,1);
plot(td,rd);
title(sprintf('\nOrder parameter magnitudes\n'));
xlim([0 Td]);
ylim([0,1]);
xlabel('time');
ylabel('r','Rotation',0);

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

% Display detrended phases PSD

subplot(2,2,4);
%plot(t,p);
plot(t,venv);
title(sprintf('\nDetrended phases\n'));
xlabel('time');
ylabel('detrended phases');
xlim([0 T]);
ylim([-1.05,+1.05]);

% Encode aggregate signal

if ~isempty(codec)
	afid = sprintf('kmutt_%03d.%s',seqno,codec);
	audiofile = fullfile(tempdir,afid);
	fprintf('\nWriting audio data to %s ...',audiofile);
	audiowrite(audiofile,y,fs);
	fprintf(' done\n\n');
	if play
		sound(y,fs);
	end
end
