# Tested on gcc (Linux 64-bit) and mingw-w64 (Windows 64-bit);
# note that MSVC (still) does not fully support the C99
# standard, so you may have problems building with it. Should
# also build okay on MacOS with gcc (not sure about clang/LLVM).

ifdef ComSpec
	OS = WIN
else
	UNAME = $(shell uname)
	ifeq ($(UNAME),Linux)
		OS = LINUX
	else
		ifeq ($(UNAME),Darwin)
			OS = MAC
		else
			OS = UNKNOWN
		endif
	endif
endif

SRC = kuramoto.c
OBJ = $(patsubst %.c,%.o,$(SRC))

MFLAGS  = -std=c99 -march=native
WFLAGS  = -Wall -Wextra -Wconversion -Winline -pedantic-errors
DFLAGS  = -D_POSIX_C_SOURCE=199309L -D_DEFAULT_SOURCE
ifeq ($(OS),WIN)
	OFLAGS    = -O3
	RM        = del /F /Q
else
	OFLAGS    = -O3 -flto
	RM        = rm -f
endif
CFLAGS = $(MFLAGS) $(OFLAGS) $(WFLAGS) $(DFLAGS)

CC  = gcc

.PHONY: all clean diag

all: $(OBJ)

$(OBJ): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJ)

diag:
	@echo "*** SRC   = " $(SRC)
