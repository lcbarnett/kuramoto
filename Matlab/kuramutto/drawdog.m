function dotfile = drawdog(nweight,eweight,gvfile,selfcon,gvprog,gvdisp,ncscale,ecscale,gscale,nwidth,nfontsize,cpendwidth,carrowsize)

% Output GraphViz 'dot' file for drawing node/edge-weighted graph, and optionally display.
%
% Requires GraphViz (https://graphviz.org/) and optionally the xdot viewer (https://github.com/jrfonseca/xdot.py)
%
% Inputs:
%
%     nweight     node weights in range [0,1], or scalar = number of nodes for unweighted nodes
%     eweight     edge weights in range [0,1] (zero for no edge)
%     gvfile      name of dot file, possibly including full path, but without extension
%     selfcon     display self-connections? (default: false)
%     gvprog      GraphViz filter ('dot, 'neato', 'fdp', etc. (see man pages; default: 'dot')
%     gvdisp      Display graph using 'xdot' (else just generate dot file; default: true)
%     nscale      node HSV scale function handle (default: blue scale)
%     escale      edge HSV scale function handle (default: red scale)

if isscalar(nweight)
	n = nweight;
	nodeweights = false;
else
	assert(isvector(nweight), 'node weights must be a scalar (number of nodes) or a vector');
	n = length(nweight);
	nodeweights = true;
end

assert(isequal(size(eweight),[n n]), 'edge weights must be a square matrix conformant with number of nodes');

if nargin < 3 || isempty(gvfile), gvfile = fullfile(tempdir,mfilename); end

if nargin < 4 || isempty(selfcon), selfcon = false; end

gvpdef = 'neato'; % also try 'dot', 'fdp'
figfmt = [];
if nargin < 5 || isempty(gvprog)
	gvprog = gvpdef;
else
	gvprog = lower(gvprog);
	k = strfind(gvprog,'/');
	if k == 1
		figfmt = gvprog(2:end);
		gvprog = gvpdef;
	elseif ~isempty(k)
		figfmt = gvprog(k+1:end);
		gvprog = gvprog(1:k-1);
	end
end

if nargin < 6 || isempty(gvdisp), gvdisp = true; end

if nodeweights
	if nargin < 7 || isempty(ncscale)
		hue = 240/360; % blue
		nHSV = [hue*ones(n,1) nweight ones(n,1)];
	elseif isscalar(nscale)
		hue = nscale;
		nHSV = [hue*ones(n,1) nweight ones(n,1)];
	elseif isa(ncscale,'function_handle')
		nHSV = ncscale(nweight);
	else
		assert(isscalar(nscale),'node color scale must be a function handle or a scalar HSV hue');
	end
else
	nHSV = [zeros(n,1) zeros(n,1) 0.85*ones(n,1)]; % light grey
end

if nargin < 8 || isempty(ecscale)
	hue = 0/360; % red
	eHSV = cat(3,hue*ones(n),eweight,ones(n));
elseif isscalar(escale)
	hue = escale;
	eHSV = cat(3,hue*ones(n),eweight,ones(n));
elseif isa(ecscale,'function_handle')
	eHSV = ecscale(eweight);
else
	assert(isscalar(nscale),'edge color scale must be a function handle or a scalar HSV hue');
end

% Default graph characteristics

if nargin <  9 || isempty(gscale),     gscale     =  2;     end
if nargin < 10 || isempty(nwidth),     nwidth     =  0.4;   end
if nargin < 11 || isempty(nfontsize),  nfontsize  =  20;    end
if nargin < 12 || isempty(cpendwidth), cpendwidth =  1.5;   end
if nargin < 13 || isempty(carrowsize), carrowsize =  1.0;   end

dotfile = [gvfile '.dot'];

df = fopen(dotfile,'wt');

fprintf(df,'digraph\n{\n');
fprintf(df,'\tgraph [fontsize = 10, scale = %g, overlap = true];\n',gscale);
fprintf(df,'\tnode [shape = circle, style = filled, fixedsize = true, width = %g, fontsize = %d];\n',nwidth,nfontsize);
fprintf(df,'\tedge [splines = true, penwidth = %g, arrowsize = %g];\n\n',cpendwidth,carrowsize);

for k = 1:n
	fprintf(df,'\t%d [fillcolor = "%6.8f,%6.8f,%6.8f"]\n',k,nHSV(k,1),nHSV(k,2),nHSV(k,3));
end
fprintf(df,'\n');

for i = 1:n
	for j=1:n
		if eweight(i,j) > 0
			if j == i
				if selfcon
					fprintf(df,'\t%d -> %d [color = "%6.8f,%6.8f,%6.8f"]\n',j,i,eHSV(i,i,1),eHSV(i,i,2),eHSV(i,i,3));
				end
			else
				fprintf(df,'\t%d -> %d [color = "%6.8f,%6.8f,%6.8f"]\n',j,i,eHSV(i,j,1),eHSV(i,j,2),eHSV(i,j,3));
			end
		end
	end
end

fprintf(df,'}\n');

fclose(df);

fprintf('drawdog: produced dot file: %s\n',dotfile);

if gvdisp
	system(sprintf('xdot -f %s %s &',gvprog,dotfile));
end

if ~isempty(figfmt)
	figfile = [gvfile '.' figfmt];
	if strcmp(figfmt,'pdf')
		epsfile = [gvfile '.eps'];
		system(sprintf('%s -Teps %s -o %s && epstopdf %s && rm %s',gvprog,dotfile,epsfile,epsfile,epsfile));
	else
		system(sprintf('%s -T%s %s -o %s',gvprog,figfmt,dotfile,figfile));
	end
	fprintf('drawdog: produced fig file: %s\n',figfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   HSV Summary                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HUE
% ---
%
% Hue is the color portion of the model, expressed as a number from 0 to 360 degrees:
%
% Red     falls between   0-60  degrees
% Yellow  falls between  61-120 degrees
% Green   falls between 121-180 degrees
% Cyan    falls between 181-240 degrees
% Blue    falls between 241-300 degrees
% Magenta falls between 301-360 degrees
%
% SATURATION
% ----------
%
% Saturation describes the amount of gray in a particular color, from 0 to 100 percent.
% Reducing this component toward zero introduces more gray and produces a faded effect.
% Sometimes, saturation appears as a range from just 0-1, where 0 is gray, and 1 is a
% primary color.
%
% VALUE (OR BRIGHTNESS)
% ---------------------
%
% Value works in conjunction with saturation and describes the brightness or intensity
% of the color, from 0-100 percent, where 0 is completely black, and 100 is the brightest
% and reveals the most color.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
