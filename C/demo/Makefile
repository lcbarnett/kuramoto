# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = main.c demo.c audio.c stulan.c chaos.c lorenz96.c kutils.c mt64.c clap.c
INC = kutils.h clap.h

BIN = kuramoto

LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto
BINDIR ?= $(PREFIX)/bin

RM ?= rm -f
CC ?= gcc

# If you don't have Gnuplot installed, omit the  -D_HAVE_GNUPLOT flag
# Deliberately NOT asking for -Winline (it is what it is)
CFLAGS = -std=c99 -march=native -O3 -flto -Wall -Wextra -Wno-unused-parameter -pedantic-errors -D_DEFAULT_SOURCE -D_HAVE_GNUPLOT


.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC) $(INC)
	$(CC) $(SRC) $(CFLAGS) -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -lkuramoto -lm -lpthread -o $(BIN)

clean:
	$(RM) $(BIN)
