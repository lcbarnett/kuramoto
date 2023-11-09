fprintf('[kuramoto startup] Initialising kuramoto Matlab interface\n');

global kuramoto_root;
kuramoto_root = fileparts(mfilename('fullpath')); % directory containing this file
addpath(kuramoto_root);
fprintf('[kuramoto startup] Added path %s and appropriate subpaths\n',kuramoto_root);

if exist('kuramoto_mex','file') ~= 3;
	fprintf(2,'[kuramoto startup] ERROR: Missing mex file ''kuramoto_mex''; please see build instructions at https://github.com/lcbarnett/kuramoto\n');
	fprintf(2,'[kuramoto startup] Bailing out\n');
	return
end

% Set up some global variables

global seqno CD_SRATE MIDDLE_C

seqno = 0;
CD_SRATE = 44100;
MIDDLE_C = 261.625565;

fprintf('[kuramoto startup] Global variables set\n');
fprintf('[kuramoto startup] Initialised (you may re-run `startup'' at any time)\n');
