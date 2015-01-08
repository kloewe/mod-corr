#-----------------------------------------------------------------------
# File    : makefile
# Contents: build pcc and tetracc programs (on Unix systems)
# Author  : Christian Borgelt
#-----------------------------------------------------------------------
SHELL       = /bin/bash
CPUINFODIR  = ../../cpuinfo/src

CC       = gcc -std=c99 -march=native
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3
# CFLAGS   = $(CFBASE) -g
INCS     = -I$(CPUINFODIR)

LD       = gcc
LDFLAGS  =
LIBS     = -lm -lrt -lpthread

HDRS     = binarize.h $(CPUINFODIR)/cpuinfo.h
OBJS     = binarize.o $(CPUINFODIR)/cpuinfo.o
PRGS     = ../bin/pcc ../bin/tetracc

#-----------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------
all:        $(PRGS)

pcc:        ../bin/pcc
	echo ""

tetracc:    ../bin/tetracc
	echo ""

../bin/pcc:         $(OBJS) pcc_main.o makefile
	$(LD) $(LDFLAGS) $(OBJS) pcc_main.o $(LIBS) -o $@

../bin/tetracc:     $(OBJS) tcc_main.o makefile
	$(LD) $(LDFLAGS) $(OBJS) tcc_main.o $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
pcc_main.o:   pcc.h $(HDRS)
pcc_main.o:   pcc.c makefile
	$(CC) $(CFLAGS) -funroll-loops $(INCS) -DPCC_MAIN -c pcc.c -o $@

tcc_main.o:   tetracc.h $(HDRS)
tcc_main.o:   tetracc.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTCC_MAIN -c tetracc.c -o $@

#-----------------------------------------------------------------------
# Main Modules
#-----------------------------------------------------------------------
pcc.o:      pcc.h $(HDRS)
pcc.o:      pcc.c makefile
	$(CC) $(CFLAGS) -funroll-loops $(INCS) -c pcc.c -o $@

tetracc.o:  tetracc.h $(HDRS)
tetracc.o:  tetracc.c makefile
	$(CC) $(CFLAGS) $(INCS) -c tetracc.c -o $@

#-----------------------------------------------------------------------
# Utility Functions
#-----------------------------------------------------------------------
binarize.o: binarize.h
binarize.o: binarize.c makefile
	$(CC) $(CFLAGS) -c binarize.c -o $@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(CPUINFODIR)/cpuinfo.o:
	cd $(CPUINFODIR); $(MAKE) cpuinfo.o

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f corr.zip corr.tar.gz; \
        zip -rq corr.zip    corr/{src,doc,bin} cpuinfo/{src,doc}; \
        tar cfz corr.tar.gz corr/{src,doc,bin} cpuinfo/{src,doc}

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) $(HOME)/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.o *~ $(PRGS)

clean:
	$(MAKE) localclean
	cd $(CPUINFODIR); $(MAKE) clean
