function cylinder_plot(t,h)

% Display phases vs. time on a cylinder

[n,N] = size(h);
assert(isvector(t) && ismatrix(h) && length(t) == n);

c = cos(h);
s = sin(h);

plot3(t,c(:,1),s(:,1));
hold on
for i = 2:N
	plot3(t,c(:,i),s(:,i));
end
hold off

xlabel('time');
