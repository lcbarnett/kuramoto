# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

SRC    = kuramoto_demo.c
LIBSRC = kuramoto.c
LIBINC = kuramoto.h
LLIB   = kuramoto

ifdef ComSpec
#	Windows
	RM = del /F /Q
	LIB = kuramoto.dll
	BIN = kuramoto_demo.exe
	BINSTDIR ?=
	LIBINSTDIR ?=
else
#	Linux, Darwin
	RM = rm -f
	LIB = libkuramoto.so
	BIN = kuramoto_demo
	BININSTDIR ?= /usr/local/bin
	LIBINSTDIR ?= /usr/local/lib
endif

CC = gcc
INSTALL = install

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_DEFAULT_SOURCE -D_GNU_SOURCE -D_GNUPLOT_HAVE_PIPE

.PHONY: all clean install uninstall

all: $(LIB) $(BIN)

$(LIB): $(LIBSRC) $(LIBINC)
	$(CC) $(LIBSRC) -shared -fPIC $(CFLAGS) -o $(LIB)

$(BIN): $(SRC) $(LIB)
	$(CC) $(SRC) $(LIBSRC) $(CFLAGS) -lm -o $(BIN)

install:
	$(INSTALL) $(BIN) $(BININSTDIR)
	$(INSTALL) $(LIB) $(LIBINSTDIR)

uninstall:
	$(RM) $(BININSTDIR)/$(BIN)
	$(RM) $(LIBINSTDIR)/$(LIB)

clean:
	$(RM) $(LIB) $(BIN)
