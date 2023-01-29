# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

SRC = kuramoto.c
INC = kuramoto.h

ifdef ComSpec
#	Windows
	RM = del /F /Q
	LIB = kuramoto.dll
	LIBDIR ?=
	INCDIR ?=
else
#	Linux, Darwin
	RM = rm -f
	LIB = libkuramoto.so
	LIBDIR ?= /usr/local/lib
	INCDIR ?= /usr/local/include/kuramoto
endif

CC = gcc
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
