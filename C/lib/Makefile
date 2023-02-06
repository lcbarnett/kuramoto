# This Makefile targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

PREFIX ?= /usr/local

SRC = kuramoto.c stulan.c
INC = kuramoto.h stulan.h

RM = rm -f
LIB = libkuramoto.so
LIBDIR ?= $(PREFIX)/lib
INCDIR ?= $(PREFIX)/include/kuramoto

CC ?= gcc
INSTALL = install -D

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Winline -pedantic-errors -D_DEFAULT_SOURCE

.PHONY: all clean install uninstall

all: $(LIB)

$(LIB): $(SRC) $(INC)
	$(CC) $(SRC) -shared -fPIC $(CFLAGS) -o $(LIB)

install:
	$(INSTALL) -m755 $(LIB) $(LIBDIR)/$(LIB)
	$(INSTALL) -m644 $(INC) $(INCDIR)/$(INC)

uninstall:
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(INCDIR)/$(INC)

clean:
	$(RM) $(LIB)
