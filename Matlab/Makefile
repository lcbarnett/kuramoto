# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

MEXT = $(shell mexext)

SRC = kuramoto_euler_mex.c kuramoto_rk4_mex.c
BIN = $(SRC:.c=.$(MEXT))

ifdef ComSpec
#	Windows: MATLAB say to use 'CFLAGS' for MinGW, but 'COMPFLAGS' seems to work; see
#   https://uk.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
	RM      = del /F /Q
	MCFLAGS = COMPFLAGS
	OBJ     = ..\C\kuramoto.o
	INCDIR  = ..\C
	LIBDIR  = ..\lib
else
#	Linux, Darwin
	RM      = rm -f
	MCFLAGS = CFLAGS
	OBJ     = ../C/kuramoto.o
	INCDIR  = ../C
	LIBDIR  = ../lib
endif

# dynamic linkage to shared library?

ifdef SOLINK
	OLINK = -L$(LIBDIR) -lkuramoto
else
	OLINK = $(OBJ)
endif

CC  = gcc

CFLAGS = -std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -I$(INCDIR) -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L

MEX = mex

MFLAGS = -O -R2018a

.PHONY: all clean diag

all: $(BIN)

clean:
	$(RM) $(BIN)

%.$(MEXT): %.c
	$(MEX) $(MFLAGS) $(MCFLAGS)="\$$CFLAGS $(CFLAGS)" $(OLINK) $<

$(BIN): $(OBJ)

diag:
	@echo "*** MEXT = " $(MEXT)
	@echo "*** SRC  = " $(SRC)
	@echo "*** OBJ  = " $(OBJ)
	@echo "*** BIN  = " $(BIN)
