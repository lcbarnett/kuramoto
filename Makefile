# Tested on gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and
# Clang/LLVM (macOS). Note that MSVC (still) does not fully
# support the C99 standard, so may have problems building.

LIBDIR = lib
ifdef ComSpec
	RM   = del /F /Q
	CP   = copy
	INC  = $(INCDIR)\kuramoto.h
	LIB  = $(LIBDIR)\libkuramoto.dll
	IDIR = $(INSTDIR)\\
else
#	Linux, Darwin
	RM   = rm -f
	CP   = cp
	INC  = $(INCDIR)/kuramoto.h
	LIB  = $(LIBDIR)/libkuramoto.so
	IDIR = $(INSTDIR)/
endif

.PHONY: all clean install uninstall diag

all:
	make -C C && make -C Matlab

clean:
	make -C C clean && make -C Matlab clean

install:
	$(CP) $(LIB) $(IDIR)

uninstall:
	$(RM) $(IDIR)$(LIB)

diag:
	@echo "*** LIBDIR   = " $(LIBDIR)
	@echo "*** LIB      = " $(LIB)
	@echo "*** INSTDIR  = " $(INSTDIR)
	@echo "*** CP       = " $(CP)
