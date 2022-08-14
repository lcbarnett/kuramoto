
N  = 10;
w  = 0.5*randn(N,1);
K  = 1/N;
h0 = pi*(2*rand(N,1)-1);
T  = 100;
dt = 0.01;

fprintf('\n');

st1 = tic;
[h1,r1,psi1,T,n] = kuramoto(N,w,K,h0,T,dt,false);
et1 = toc(st1);
fprintf('Euler method : %g seconds\n',et1);

st2 = tic;
[h2,r2,psi2,T,n] = kuramoto(N,w,K,h0,T,dt,true);
et2 = toc(st2);
fprintf('Runge-Kutta  : %g seconds\n',et2);

mad = max(abs(r1-r2));
fprintf('\nMaximum absolute difference = %g\n\n',mad);

t = linspace(0,T,n);

figure(1); clf
plot(t',[r1;r2]');
legend({'Euler','RK4'})
xlabel('time')
ylabel('order parameter')
title('Kuramoto system')
