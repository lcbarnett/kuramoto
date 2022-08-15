# Tested on gcc (Linux 64-bit) and mingw-w64 (Windows 64-bit);
# note that MSVC (still) does not fully support the C99
# standard, so you may have problems building with it. Should
# also build okay on MacOS with gcc (not sure about clang/LLVM).

ifdef ComSpec
	OS = WIN
else
	UNAME = $(shell uname)
	ifeq ($(UNAME),Linux)
		OS = LINUX
	else
		ifeq ($(UNAME),Darwin)
			OS = MAC
		else
			OS = UNKNOWN
		endif
	endif
endif

MEXEXT =  $(shell mexext)

SRC = kuramoto_euler_mex.c kuramoto_rk4_mex.c
BIN = $(patsubst %.c,%.$(MEXEXT),$(SRC))
OBJ = ../C/kuramoto.o

MFLAGS = -std=c99 -march=native
WFLAGS = -Wall -Wextra -Wconversion -Winline -pedantic-errors
IFLAGS = -I../C
DFLAGS = -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L
OFLAG  = -O3
ifeq ($(OS),WIN)
#	MATLAB () say to use 'CFLAGS' for MinGW, but 'COMPFLAGS' seems to work
#   See: https://uk.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
	MCFLAGS = COMPFLAGS
	RM      = del /F /Q
	RMBIN   = $(subst /,\,$(BIN))
else
	MCFLAGS = CFLAGS
	RM      = rm -f
	RMBIN   = $(BIN)
endif
GCCFLAGS = $(MFLAGS) $(OFLAGS) $(WFLAGS) $(IFLAGS) $(DFLAGS)
MEXFLAGS = -O -R2018a

CC  = gcc
MEX = mex

.PHONY: all clean diag

all: $(BIN)

clean:
	$(RM) $(RMBIN)

diag:
	@echo "*** SRC   = " $(SRC)
	@echo "*** BIN   = " $(BIN)
	@echo "*** MEXEXT = " $(MEXEXT)

%.$(MEXEXT): %.c
	$(MEX) $(MEXFLAGS) $(MCFLAGS)="\$$CFLAGS $(GCCFLAGS)" $(OBJ) $<
