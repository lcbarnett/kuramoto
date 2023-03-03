# Makefiles targets POSIX environments, probably using GCC,
# but should also work with Clang/LLVM on Nac OS. It may also
# work on Windows, say with mingw-w64/MSYS2 (note that MSVC does
# not fully support the C99 standard, so may be problematic).

.PHONY: all clean install uninstall

all:
	make -C lib && make -C demo

install:
	make -C lib install && make -C demo install

uninstall:
	make -C lib uninstall && make -C demo uninstall

clean:
	make -C lib clean && make -C demo clean
