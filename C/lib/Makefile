# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = kuramoto.c kuramoto_old.c
INC = kuramoto.h kuramoto_old.h

LIB = libkuramoto.so
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto

RM ?= rm -f
CC ?= gcc

CFLAGS = -std=c99 -march=native -O3 -flto -Wall -Wextra -Winline -Wno-unused-parameter -pedantic-errors -D_DEFAULT_SOURCE

.PHONY: all clean install uninstall

all: $(LIB)

$(LIB): $(SRC) $(INC)
	$(CC) $(SRC) -shared -fPIC $(CFLAGS) -o $(LIB)

install:
	install -m755 $(LIB) $(LIBDIR)
	install -D -m644 kuramoto.h $(INCDIR)/kuramoto.h
	install -D -m644 kuramoto_old.h $(INCDIR)/kuramoto_old.h

uninstall:
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(INCDIR)/kuramoto.h
	$(RM) $(INCDIR)/kuramoto_old.h

clean:
	$(RM) $(LIB)
