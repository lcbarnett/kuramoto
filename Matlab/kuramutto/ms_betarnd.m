function x = ms_betarnd(varargin)

mu    = varargin{1};
sigma = varargin{2};

g = ((mu*(1-mu))/sigma^2)-1;
a = g*mu;
b = g*(1-mu);

x = betarnd(a,b,varargin{3:end});
