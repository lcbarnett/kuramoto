# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

SRC = kuramoto.c
OBJ = $(SRC:.c=.o)

CFLAGS = -std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L

ifdef ComSpec
#	Windows
	LIBDIR = ..\lib
	LIB = $(LIBDIR)\libkuramoto.dll
	RM = del /F /Q
else
#	Linux, Darwin
	LIBDIR = ../lib
	LIB = $(LIBDIR)/libkuramoto.so
	RM = rm -f
endif

CC  = gcc

.PHONY: all clean diag

all: $(LIB)

$(LIB): $(OBJ)
	$(CC) -shared -o $@ $^


$(OBJ): %.o: %.c
	$(CC) -c -fPIC $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJ) $(LIB)

diag:
	@echo "*** SRC = " $(SRC)
	@echo "*** OBJ = " $(OBJ)

kuramoto.o: kuramoto.h
