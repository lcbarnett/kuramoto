% Kuramutto connectivity spec

% Default parameters (override on command line - see 'defvar.m')

defvar('kmean', 100     ); % oscillator coupling constants mean      (frequency multiplier - dimensionless)
defvar('ksdev', kmean/2 ); % oscillator coupling constants std. dev. (frequency multiplier - dimensionless)
defvar('kinhp', 0.1     ); % oscillator coupling probability of inhibitory connection
defvar('kseed', []      ); % oscillator couplinghttps://www.mathworks.com/help/releases/R2023a/matlab/ref/audiowrite.html?doclanguage=en&nocookie=true&prodfilter=ML%20SL%205G%20AE%20AT%20AA%20AU%20DR%20AS%20BI%20BL%20C2%20CM%20VP%20CT%20CF%20DB%20DF%20DD%20DH%20noscnosc%20HS%20DS%20ET%20EC%20FH%20IT%20FI%20PO%20FL%20GD%20GC%20HD%20ES%20IA%20IP%20OT%20IC%20LP%20LS%20MG%20ME%20CO%20MJ%20MR%20TE%20DX%20AM%20MP%20MT%20noscV%20OP%20DM%20PD%20AR%20PW%20PM%20RA%20RL%20RQ%20RB%20RP%20RF%20RK%20RO%20RC%20RR%20SI%20TF%20SX%20SQ%20SG%20SB%20SE%20SS%20BT%20LD%20PS%20SH%20MS%20VR%20VV%20CI%20RT%20SK%20SD%20CV%20SO%20DV%20PL%20XP%20SR%20SZ%20HW%20SF%20ST%20SM%20ZC%20ID%20TA%20UV%20VE%20Vnosc%20VT%20WA%20LH%20WB%20WL&docviewer=helpbrowser&docrelease=R2023a&s_cid=pl_webdoc&loginurl=https%3A%2F%2F127.0.0.1%3A31515%2Ftoolbox%2Fmatlab%2Flogin%2Fweb%2Findex.html%3Fsnc%3DIAnoscGSD%26external%3Dtrue%26channel%3D__mlfpmc__&searchsource=mw&snc=6TZA19&container=jshelpbrowser#d124e12665 constants random seed (empty for no seeding)

if ~isempty(kseed), rstate = rng(kseed); end

K = ms_lognrnd(kmean,ksdev,nosc,nosc);

kinhib = rand(nosc)<kinhp;
K(kinhib) = -K(kinhib);

K = C.*K;

if ~isempty(kseed), rng(rstate); end

fprintf('Oscillator coupling constants =\n\n');
disp(K)
