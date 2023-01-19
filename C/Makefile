# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

SRC = kuramoto.c
OBJ = $(SRC:.c=.o)

CFLAGS = -std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_DEFAULT_SOURCE

ifdef ComSpec
#	Windows
	RM = del /F /Q
	LIBDIR = ..\lib
	LIB = $(LIBDIR)\libkuramoto.dll
else
#	Linux, Darwin
	RM = rm -f
	LIBDIR = ../lib
	LIB = $(LIBDIR)/libkuramoto.so
endif

CC  = gcc

.PHONY: all clean diag

all: $(LIB)

# Build shared library

$(LIB): $(OBJ)
	$(CC) -shared -o $@ $^

$(OBJ): %.o: %.c
	$(CC) -c -fPIC $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJ) $(LIB)

kuramoto.o: kuramoto.h

diag:
	@echo "*** SRC = " $(SRC)
	@echo "*** OBJ = " $(OBJ)
