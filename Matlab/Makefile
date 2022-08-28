# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

MEXT = $(shell mexext)

SRC = kuramoto_euler_mex.c kuramoto_rk4_mex.c kuramoto_noisy.c
BIN = $(patsubst %.c,%.$(MEXT),$(SRC))
OBJ = ../C/kuramoto.o

CFLAGS = -std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -I../C -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L

MFLAGS = -O -R2018a

ifdef ComSpec
#	Windows: MATLAB say to use 'CFLAGS' for MinGW, but 'COMPFLAGS' seems to work; see
#   https://uk.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
	MCFLAGS = COMPFLAGS
	RM      = del /F /Q
	RMBIN   = $(subst /,\,$(BIN))
else
#	Linux, Darwin
	MCFLAGS = CFLAGS
	RM      = rm -f
	RMBIN   = $(BIN)
endif

CC  = gcc
MEX = mex

.PHONY: all clean diag

all: $(BIN)

clean:
	$(RM) $(RMBIN)

diag:
	@echo "*** SRC  = " $(SRC)
	@echo "*** BIN  = " $(BIN)
	@echo "*** OBJ  = " $(OBJ)
	@echo "*** MEXT = " $(MEXT)

%.$(MEXT): %.c
	$(MEX) $(MFLAGS) $(MCFLAGS)="\$$CFLAGS $(CFLAGS)" $(OBJ) $<
