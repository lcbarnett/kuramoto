# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

SRC = kuramoto.c
OBJ = $(patsubst %.c,%.o,$(SRC))

CFLAGS = -std=c99 -march=native -O3 -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L

ifdef ComSpec
#	Windows
	RM = del /F /Q
else
#	Linux, Darwin
	RM = rm -f
endif

CC  = gcc

.PHONY: all clean diag

all: $(OBJ)

$(OBJ): %.o: %.c
	$(CC) -c -fPIC $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJ)

diag:
	@echo "*** OBJ = " $(SRC)
	@echo "*** OBJ = " $(SRC)

kuramoto.o: kuramoto.h
