function x = ms_lognrnd(varargin)

mu    = varargin{1};
sigma = varargin{2};

u = log(1+sigma^2/mu^2);

x = lognrnd(log(mu)-u/2,sqrt(u),varargin{3:end});
