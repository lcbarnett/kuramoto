
% Default parameters (override on command line - see 'defvar.h')

defvar('N',     20      ); % number of oscillators
defvar('wmean', 0       ); % oscillator frequencies mean
defvar('wsdev', pi/7    ); % oscillator frequencies std. dev.
defvar('wseed', []      ); % oscillator frequencies random seed (empty for no seeding)
defvar('Kmean', 0.8/N   ); % oscillator coupling constants mean
defvar('Ksdev', Kmean/6 ); % oscillator coupling constants std. dev.
defvar('Kseed', []      ); % oscillator coupling constants random seed (empty for no seeding)
defvar('a',     0       ); % oscillator phase lag constant
defvar('hseed', []      ); % oscillator initial phases random seed (empty for no seeding)
defvar('T',     200     ); % simulation time
defvar('dt',    0.01    ); % integration time increment
defvar('nmean', 0.05    ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev', nmean/3 ); % oscillator input noise magnitude std. dev.
defvar('nseed', []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed', []      ); % oscillator input noise random seed (empty for no seeding)

% Random Kuramoto parameters

if ~isempty(wseed), rstate = rng(wseed); end
w = wmean + wsdev*randn(N,1); % oscillator frequencies normally distributed with mean wmean and std. dev wsdev
if ~isempty(wseed), rng(rstate); end

if ~isempty(Kseed), rstate = rng(Kseed); end
K = Kmean + Ksdev*randn(N); % coupling constants normally distributed with mean Kmean and std. dev. Ksdev
if ~isempty(Kseed), rng(rstate); end

if ~isempty(hseed), rstate = rng(hseed); end
h0 = pi*(2*rand(N,1)-1);      % initial phases uniform on [-pi,pi]
if ~isempty(hseed), rng(rstate); end

n = round(T/dt);
assert(n > 0,'Simulation time too short, or time increment too large!');
T = n*dt; % adjusted simulation time

if nmean > 0
	if ~isempty(nseed), rstate = rng(nseed); end
	nmag = gamrnd(nmean^2/nsdev^2,nsdev^2/nmean,N,1); % per-oscillator noise magnitudes drawn from Gamma distribution
	if ~isempty(nseed), rng(rstate); end

	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = nmag.*randn(N,n); % uncorrelated Gaussian white noise inputs
	if ~isempty(Iseed), rng(rstate); end
else
	I = []; % no input
end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('\n');

st1 = tic;
[h1,r1,psi1] = kuramoto(N,w,K,a,h0,n,dt,'Euler',I);
et1 = toc(st1);
fprintf('Euler method : %g seconds\n',et1);

st2 = tic;
[h2,r2,psi2] = kuramoto(N,w,K,a,h0,n,dt,'RK4');
et2 = toc(st2);
fprintf('RK4 method   : %g seconds\n',et2);

mad = max(abs(r1-r2));
fprintf('\nMaximum absolute difference = %g\n\n',mad);

t = linspace(0,T,n);

% Display order parameter magnitude

figure(1); clf;
plot(t',[r1;r2]');
legend({'Euler','RK4'});
xlabel('time');
ylabel('order parameter magnitude (r)');
title(sprintf('\nKuramoto system: N = %d\n',N));

% Display oscillator phases (RK4) on cylinder

cosh1 = cos(h1);
sinh1 = sin(h1);
cosh2 = cos(h2);
sinh2 = sin(h2);

figure(2); clf;
sgtitle(sprintf('\nKuramoto system: N = %d\n',N));
subplot(2,1,1);
plot3(t',cosh1(1,:)',sinh1(1,:)');
xlabel('time');
ylabel('x');
zlabel('y');
hold on
for i = 2:N
	plot3(t',cosh1(i,:)',sinh1(i,:)');
end
hold off
title(sprintf('Euler\n'));
subplot(2,1,2);
plot3(t',cosh2(1,:)',sinh2(1,:)');
xlabel('time');
ylabel('x');
zlabel('y');
hold on
for i = 2:N
	plot3(t',cosh2(i,:)',sinh2(i,:)');
end
hold off
title(sprintf('RK4\n'));

% Display order parameter (animation)

x1 = r1.*cos(psi1);
y1 = r1.*sin(psi1);
x2 = r2.*cos(psi2);
y2 = r2.*sin(psi2);

figure(3); clf;
rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
daspect([1,1,1]);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
set(gca,'XTick',[])
set(gca,'YTick',[])
box on
title(sprintf('\nOrder parameter (z)\n'));
ln1 = line([0;0],[0;0],'LineWidth',1,'Color','b');
ln2 = line([0;0],[0;0],'LineWidth',1,'Color','r');
ts = xlabel('t = 0');
legend({'Euler','RK4'});
for k = 1:n
	ln1.XData = [0;x1(k)];
	ln1.YData = [0;y1(k)];
	ln2.XData = [0;x2(k)];
	ln2.YData = [0;y2(k)];
	ts.String = sprintf('t = %.0f',k*dt);
	drawnow limitrate nocallbacks
end
