# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = kuramoto.c stulan.c ode.c
INC = kuramoto.h stulan.h ode.h

LIB = libkuramoto.so
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto

RM ?= rm -f
CC ?= gcc

# Deliberately NOT asking for -Winline (it is what it is)
CFLAGS = -std=c99 -march=native -O3 -flto -Wall -Wextra -Wno-unused-parameter -pedantic-errors -D_DEFAULT_SOURCE

.PHONY: all clean install uninstall

all: $(LIB)

$(LIB): $(SRC) $(INC)
	$(CC) $(SRC) -shared -fPIC $(CFLAGS) -o $(LIB)

install:
	install -m755 $(LIB) $(LIBDIR)
	install -D -m644 kuramoto.h $(INCDIR)/kuramoto.h
	install -D -m644 stulan.h $(INCDIR)/stulan.h
	install -D -m644 ode.h $(INCDIR)/ode.h

uninstall:
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(INCDIR)/kuramoto.h
	$(RM) $(INCDIR)/stulan.h
	$(RM) $(INCDIR)/ode.h

clean:
	$(RM) $(LIB)
