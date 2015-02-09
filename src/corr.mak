#-----------------------------------------------------------------------
# File    : corr.mak
# Contents: build pcc and tetracc programs (on Windows systems)
# Author  : Christian Borgelt
#-----------------------------------------------------------------------
CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS
ARCH     = /D __SSE2__ /D __SSSE3__ # /D __SSE4_1__ /D __POPCNT__
CFLAGS   = /nologo /W3 /O2 /GS- $(DEFS) $(ARCH) /c $(ADDFLAGS)
INCS     =

LD       = link.exe
LDFLAGS  = /nologo /subsystem:console /incremental:no
LIBS     = winmm.lib

HDRS     = binarize.h   cpuinfo.h
OBJS     = binarize.obj cpuinfo.obj
PRGS     = pcc.exe tetracc.exe

#-----------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------
all:          $(PRGS)

pcc.exe:      $(OBJS) pcc.obj pcc.mak
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) pcc.obj     /out:$@

tetracc.exe:  $(OBJS) tetracc.obj pcc.mak
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) tetracc.obj /out:$@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
pcc.obj:      pcc.h $(HDRS)
pcc.obj:      pcc.c pcc.mak
	$(CC) $(CFLAGS) $(INCS) /D PCC_MAIN pcc.c     /Fo$@

tetracc.obj:  tetracc.h $(HDRS)
tetracc.obj:  tetracc.c tetracc.h pcc.mak
	$(CC) $(CFLAGS) $(INCS) /D TCC_MAIN tetracc.c /Fo$@

#-----------------------------------------------------------------------
# Utility Functions
#-----------------------------------------------------------------------
binarize.obj: $(HDRS)
binarize.obj: binarize.c makefile
	$(CC) $(CFLAGS) $(INCS) binarize.c /Fo$@

cpuinfo.obj:  $(HDRS)
cpuinfo.obj:  cpuinfo.c makefile
	$(CC) $(CFLAGS) $(INCS) cpuinfo.c /Fo$@

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	-@copy pcc.exe     ..\..\..\bin
	-@copy tetracc.exe ..\..\..\bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	-@erase /Q *~ *.obj *.idb *.pch $(PRGS)

clean:
	$(MAKE) /f pcc.mak localclean
