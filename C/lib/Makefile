# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = kuramoto.c stulan.c
INC = kuramoto.h stulan.h

LIB = libkuramoto.so
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto

CC ?= gcc

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Winline -pedantic-errors -D_DEFAULT_SOURCE

.PHONY: all clean install uninstall

all: $(LIB)

$(LIB): $(SRC) $(INC)
	$(CC) $(SRC) -shared -fPIC $(CFLAGS) -o $(LIB)

install:
	install -m755 $(LIB) $(LIBDIR)
	install -D -m644 kuramoto.h $(INCDIR)/kuramoto.h
	install -D -m644 stulan.h   $(INCDIR)/stulan.h

uninstall:
	rm $(LIBDIR)/$(LIB)
	rm $(INCDIR)/kuramoto.h
	rm $(INCDIR)/stulan.h

clean:
	rm $(LIB)
