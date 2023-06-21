fprintf('[kuramoto startup] Initialising kuramoto Matlab interface\n');

% Set up some global variables

global seqno CD_SRATE MIDDLE_C

seqno = 0;
CD_SRATE = 44100;
MIDDLE_C = 261.625565;

fprintf('[kuramoto startup] Global variables set\n');
fprintf('[kuramoto startup] Initialised (you may re-run `startup'' at any time)\n');
