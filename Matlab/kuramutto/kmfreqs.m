% Kuramutto frequency spec

if ~isempty(wheadseed), rstate = rng(wheadseed); end
whead = wheadmin+(wheadmax-wheadmin)*rand(nhead,1);
if ~isempty(wheadseed), rng(rstate); end
fprintf('\nHead natural frequencies (Hz) = '); disp(whead')

if ~isempty(wbodiseed), rstate = rng(wbodiseed); end
wbodi = wbodimin+(wbodimax-wbodimin)*rand(nbodi,1);
if ~isempty(wbodiseed), rng(rstate); end
fprintf('\nHead natural frequencies (Hz) = '); disp(wbodi')

if ~isempty(wflegseed), rstate = rng(wflegseed); end
wfleg = wflegmin+(wflegmax-wflegmin)*rand(nfleg,1);
if ~isempty(wflegseed), rng(rstate); end
fprintf('\nHead natural frequencies (Hz) = '); disp(wfleg')

if ~isempty(whlegseed), rstate = rng(whlegseed); end
whleg = whlegmin+(whlegmax-whlegmin)*rand(nhleg,1);
if ~isempty(whlegseed), rng(rstate); end
fprintf('\nHead natural frequencies (Hz) = '); disp(whleg')

if ~isempty(wtailseed), rstate = rng(wtailseed); end
wtail = wtailmin+(wtailmax-wtailmin)*rand(ntail,1);
if ~isempty(wtailseed), rng(rstate); end
fprintf('\nHead natural frequencies (Hz) = '); disp(wtail')

w = [whead; wbodi; wfleg; whleg; wtail];
