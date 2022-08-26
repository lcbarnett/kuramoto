function clock_plot(x,y,ptitle,legend,cols)

% In progress

[nvars,nphases] = size(x);

rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
daspect([1,1,1]);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
set(gca,'XTick',[])
set(gca,'YTick',[])
box on
title(ptitle);
ln = cell(nvars,1);
for i = 1:nvars
	ln{i} = line([0;0],[0;0],'LineWidth',1,'Color',col{i});
end
ts = xlabel('t = 0');
legend(plegend);
for k = 1:nphases
	for i = 1:nvars
		ln{i}.XData = [0;x{i}(k)];
		ln{i}.YData = [0;y{i}(k)];
	end
	ts.String = sprintf('t = %.0f',k*dt);
	drawnow limitrate nocallbacks
end
