% Uncorrelated Gaussian white inoise input to tail, magnitudes log-normally distributed

% Default parameters (override on command line - see 'defvar.m')

defvar('mmean', 0.01    ); % oscillator input inoise magnitude mean (zero for no inoise)
defvar('msdev', mmean/5 ); % oscillator input inoise magnitude std. dev.
defvar('mseed', []      ); % oscillator input inoise magnitude random seed (empty for no seeding)
defvar('nseed', []      ); % oscillator input inoise random seed (empty for no seeding)

if mmean > 0 % with input inoise
	if ~isempty(mseed), rstate = rng(mseed); end
	nln = log(1+msdev^2/mmean^2);
	nmag = (lognrnd(log(mmean)-nln/2,sqrt(nln),1,ntail));
	if ~isempty(mseed), rng(rstate); end

	if ~isempty(nseed), rstate = rng(nseed); end
	inoise = zeros(ndt,nosc);
	inoise(:,itail) = nmag.*randn(ndt,ntail);
	if ~isempty(nseed), rng(rstate); end
else
	inoise = []; % no inoise
end
