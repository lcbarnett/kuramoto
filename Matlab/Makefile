# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.
# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

MEXT = $(shell mexext)

SRC = kuramoto_mex.c
BIN = kuramoto.$(MEXT)

RM      = rm -f
MCFLAGS = CFLAGS
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto

CC ?= gcc

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Winline -pedantic-errors -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -D_DEFAULT_SOURCE

MEX = mex

MFLAGS = -O -R2018a

.PHONY: all clean diag

all: $(BIN)

clean:
	$(RM) $(BIN)

$(BIN): $(SRC)
	$(MEX) $(MFLAGS) $(MCFLAGS)="\$$CFLAGS $(CFLAGS)" -lkuramoto $<

diag:
	@echo "*** MEXT = " $(MEXT)
	@echo "*** SRC  = " $(SRC)
	@echo "*** BIN  = " $(BIN)
