% Uncorrelated Gaussian white inoise input to legs, magnitudes log-normally distributed

% Default parameters (override on command line - see 'defvar.m')

defvar('mmean', 0.01    ); % oscillator input noise magnitude mean (zero for no inoise)
defvar('msdev', mmean/5 ); % oscillator input noise magnitude std. dev.
defvar('mseed', []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('nseed', []      ); % oscillator input noise random seed (empty for no seeding)

ilegs = [ifleg ihleg];
nlegs = nfleg+nhleg;

if mmean > 0 % with input inoise
	if ~isempty(mseed), rstate = rng(mseed); end
	nln = log(1+msdev^2/mmean^2);
	nmag = (lognrnd(log(mmean)-nln/2,sqrt(nln),1,nlegs));
	if ~isempty(mseed), rng(rstate); end

	if ~isempty(nseed), rstate = rng(nseed); end
	inoise = zeros(ndt,nosc);
	inoise(:,ilegs) = nmag.*randn(ndt,nlegs);
	if ~isempty(nseed), rng(rstate); end
else
	inoise = []; % no inoise
end
