function y = ldetrend(x,t)

% Subtract linear fit of t to x from x

assert(length(t) == size(x,1));
tm = mean(t);
tv = mean(t.^2);
xm = mean(x);
tx = mean(t.*x);
y = x - (tv*xm-tm*tx+(tx-tm*xm).*t)/(tv-tm^2);
