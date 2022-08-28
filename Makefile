# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

.PHONY: all clean

all:
	make -C C && make -C Matlab

clean:
	make -C C clean && make -C Matlab clean
