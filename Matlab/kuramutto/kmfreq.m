% Kuramutto frequency spec

global MIDDLE_C

% Default parameters (override on command line - see 'defvar.m')

defvar('wheadseed', []         );
defvar('wheadmin',  50         );
defvar('wheadmax',  2*MIDDLE_C );
defvar('wbodiseed', []         );
defvar('wbodimin',  20         );
defvar('wbodimax',  3*MIDDLE_C );
defvar('wflegseed', []         );
defvar('wflegmin',  20         );
defvar('wflegmax',  100        );
defvar('whlegseed', []         );
defvar('whlegmin',  1          );
defvar('whlegmax',  10         );
defvar('wtailseed', []         );
defvar('wtailmin',  0          );
defvar('wtailmax',  5          );

if ~isempty(wheadseed), rstate = rng(wheadseed); end
whead = wheadmin+(wheadmax-wheadmin)*rand(nhead,1);
if ~isempty(wheadseed), rng(rstate); end
fprintf('Head natural frequencies (Hz) = '); disp(whead')

if ~isempty(wbodiseed), rstate = rng(wbodiseed); end
wbodi = wbodimin+(wbodimax-wbodimin)*rand(nbodi,1);
if ~isempty(wbodiseed), rng(rstate); end
fprintf('Body natural frequencies (Hz) = '); disp(wbodi')

if ~isempty(wflegseed), rstate = rng(wflegseed); end
wfleg = wflegmin+(wflegmax-wflegmin)*rand(nfleg,1);
if ~isempty(wflegseed), rng(rstate); end
fprintf('Front leg natural frequencies (Hz) = '); disp(wfleg')

if ~isempty(whlegseed), rstate = rng(whlegseed); end
whleg = whlegmin+(whlegmax-whlegmin)*rand(nhleg,1);
if ~isempty(whlegseed), rng(rstate); end
fprintf('Hind leg natural frequencies (Hz) = '); disp(whleg')

if ~isempty(wtailseed), rstate = rng(wtailseed); end
wtail = wtailmin+(wtailmax-wtailmin)*rand(ntail,1);
if ~isempty(wtailseed), rng(rstate); end
fprintf('Tail natural frequencies (Hz) = '); disp(wtail')

w = [whead; wbodi; wfleg; whleg; wtail];
