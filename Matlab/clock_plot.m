function clock_plot(t,r,psi,leg,col)

% Animate (complex) phases vs time

[n,m] = size(r);
assert(isvector(t) && ismatrix(r) && ismatrix(psi) && length(t) == n && size(psi,1) == n && size(psi,2) == m);
assert(isvector(leg) && length(leg) == m && isvector(col) && length(col) == m);

x = r.*cos(psi);
y = r.*sin(psi);

rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
daspect([1,1,1]);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
set(gca,'XTick',[])
set(gca,'YTick',[])
box on
ln = cell(m,1);
for i = 1:m
	ln{i} = line([0;0],[0;0],'LineWidth',1,'Color',col{i});
end
legend(leg);
ts = xlabel('t = 0');
for k = 1:n
	for i = 1:m
		ln{i}.XData = [0;x(k,i)];
		ln{i}.YData = [0;y(k,i)];
	end
	ts.String = sprintf('t = %.0f',t(k));
	drawnow limitrate nocallbacks
end
