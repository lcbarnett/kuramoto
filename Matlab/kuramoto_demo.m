
% Default parameters (override on command line - see 'defvar.h')

defvar('N',     20     ); % number of oscillators
defvar('wmean', 0      ); % oscillator frequencies mean
defvar('wsdev', pi/7   ); % oscillator frequencies std. dev.
defvar('wseed', []     ); % oscillator frequencies random seed (empty for no seeding)
defvar('Kmean', 0.8/N  ); % oscillator coupling constants mean
defvar('Ksdev', 0.1/N  ); % oscillator coupling constants std. dev.
defvar('Kseed', []     ); % oscillator coupling constants random seed (empty for no seeding)
defvar('a',     0      ); % oscillator phase lag constant
defvar('Vmean', 0.04   ); % oscillator input noise mean (zero for no noise)
defvar('Vsdev', 0.01   ); % oscillator input noise std. dev.
defvar('Vseed', []     ); % oscillator input noise random seed (empty for no seeding)
defvar('hseed', []     ); % oscillator initial phases random seed (empty for no seeding)
defvar('T',     200    ); % simulation time
defvar('dt',    0.01   ); % integration time increment

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

if Vmean > 0
	mode = 'Noisy';
	if ~isempty(Vseed), rstate = rng(Vseed); end
	V = gamrnd(Vmean^2/Vsdev^2,Vsdev^2/Vmean,N,1); % Gamma noise variance
	if ~isempty(Vseed), rng(rstate); end
else
	mode = 'Euler';
	V = mode;
end

% Run Kuramoto Euler and Rung-Kutta simulations with specified parameters

fprintf('\n');

st1 = tic;
[h1,r1,psi1,T,n] = kuramoto(N,w,K,a,h0,T,dt,V);
et1 = toc(st1);
fprintf('%s method : %g seconds\n',mode,et1);

st2 = tic;
[h2,r2,psi2,T,n] = kuramoto(N,w,K,a,h0,T,dt,'RK4');
et2 = toc(st2);
fprintf('Runge-Kutta  : %g seconds\n',et2);

mad = max(abs(r1-r2));
fprintf('\nMaximum absolute difference = %g\n\n',mad);

t = linspace(0,T,n);

% Display order parameter magnitude

figure(1); clf;
plot(t',[r1;r2]');
legend({mode,'RK4'});
xlabel('time');
ylabel('order parameter magnitude (r)');
title(sprintf('\nKuramoto system: N = %d\n',N));

% Display oscillator phases (RK4) on cylinder

cosh2 = cos(h2);
sinh2 = sin(h2);

figure(2); clf;
plot3(t',cosh2(1,:)',sinh2(1,:)');
xlabel('time');
ylabel('x');
zlabel('y');
hold on
for i = 2:N
	plot3(t',cosh2(i,:)',sinh2(i,:)');
end
hold off
title(sprintf('\nKuramoto system: N = %d\n',N));

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
legend({mode,'RK4'});
for k = 1:n
	ln1.XData = [0;x1(k)];
	ln1.YData = [0;y1(k)];
	ln2.XData = [0;x2(k)];
	ln2.YData = [0;y2(k)];
	ts.String = sprintf('t = %.0f',k*dt);
	drawnow limitrate nocallbacks
end
