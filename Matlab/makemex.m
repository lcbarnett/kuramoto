function makemex

% Build mex files - make sure you have run 'make' in the C subdirectory first!
%
% Seems to work on Linux, MacOS (with gcc or Clang) and Windows with MinGW

if ispc % Windows
	cfstr = 'COMPFLAGS'; % or maybe just 'CFLAGS'? Seems to depend on Matlab version
else    % Linux (gcc), MacOS (gcc or Clang)
	cfstr = 'CFLAGS';
end

CFLAGS = '-std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -I../C -D_POSIX_C_SOURCE=199309L -D_DEFAULT_SOURCE';

mexinvoc = ['mex -O -R2018a ' cfstr '="\$' cfstr ' ' CFLAGS '" ../C/kuramoto.o'];

eval([mexinvoc ' kuramoto_euler_mex.c']);
eval([mexinvoc ' kuramoto_rk4_mex.c']);
