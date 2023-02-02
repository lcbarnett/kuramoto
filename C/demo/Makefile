# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). We basically assume Posix compliance. Note
# that MSVC (still) does not fully support the C99 standard, so
# may have problems building.
#
# NOTE: you must build and install the kuramoto dynamic shared
# library first! See ../lib

SRC = main.c demo.c audio.c kutils.c clap.c
INC = kutils.h clap.h

ifdef ComSpec
#	Windows
	RM = del /F /Q
	BIN = kuramoto.exe
	# I don't know standard locations for dev headers and dlls on Windows
	LIBDIR ?=
	INCDIR ?=
else
#	Linux, Darwin
	RM = rm -f
	BIN = kuramoto
	LIBDIR ?= /usr/local/lib
	INCDIR ?= /usr/local/include/kuramoto
endif

CC = gcc
INSTALL = install

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Winline -pedantic-errors -D_DEFAULT_SOURCE -D_GNU_SOURCE -D_HAVE_GNUPLOT

.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC) $(INC)
	$(CC) $(SRC) $(LIBSRC) $(CFLAGS) -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -lkuramoto -lm -o $(BIN)

clean:
	$(RM) $(BIN)
