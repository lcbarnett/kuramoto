# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = main.c demo.c audio.c kutils.c clap.c
INC = kutils.h clap.h

RM = rm -f
BIN = kuramoto

LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto
BINDIR ?= $(PREFIX)/bin

CC ?= gcc

# If you don't have Gnuplot installed, omit the  -D_HAVE_GNUPLOT flag
CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Winline -pedantic-errors -D_DEFAULT_SOURCE -D_HAVE_GNUPLOT

.PHONY: all clean install uninstall

all: $(BIN)

$(BIN): $(SRC) $(INC)
	$(CC) $(SRC) $(CFLAGS) -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -lkuramoto -lm -o $(BIN)

install:
	install -D -m755 $(BIN) $(BINDIR)/$(BIN)

uninstall:
	$(RM) $(BINDIR)/$(BIN)

clean:
	$(RM) $(BIN)
