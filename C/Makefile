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
	LIBDIR = ..\lib
	LIB = $(LIBDIR)\kuramoto.dll
	BIN = kuramoto_demo.exe
else
#	Linux, Darwin
	RM = rm -f
	LIBDIR = ../lib
	LIB = $(LIBDIR)/libkuramoto.so
	BIN = kuramoto_demo
endif

CC = gcc

CFLAGS = -std=c99 -march=native -O2 -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_DEFAULT_SOURCE

.PHONY: all clean

all: $(LIB) $(BIN)

$(LIB): $(LIBSRC) $(LIBINC)
	$(CC) $(LIBSRC) -shared -fPIC $(CFLAGS) -o $(LIB)

$(BIN): $(SRC) $(LIB)
	$(CC) $(SRC) $(CFLAGS) -L$(LIBDIR) -Wl,-rpath=$(LIBDIR) -l$(LLIB) -lm -o $(BIN)

clean:
	$(RM) $(LIB) $(BIN)
