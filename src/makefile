#-----------------------------------------------------------------------
# File    : makefile
# Contents: build pcc and tetracc programs (on Unix systems)
# Author  : Christian Borgelt
#-----------------------------------------------------------------------
SHELL       = /bin/bash
CPUINFODIR  = ../../cpuinfo/src

CC       = gcc -std=c99 -march=native
# CC       = gcc -std=c99 -mavx
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

../bin/pcc:         $(OBJS) pcc.o makefile
	$(LD) $(LDFLAGS) $(OBJS) pcc.o     $(LIBS) -o $@

../bin/tetracc:     $(OBJS) tetracc.o makefile
	$(LD) $(LDFLAGS) $(OBJS) tetracc.o $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
pcc.o:      pcc.h $(HDRS)
pcc.o:      pcc.c makefile
	$(CC) $(CFLAGS) -funroll-loops $(INCS) -DPCC_MAIN -c pcc.c -o $@

tetracc.o:  tetracc.h $(HDRS)
tetracc.o:  tetracc.c tetracc.h makefile
	$(CC) $(CFLAGS) $(INCS) -DTCC_MAIN -c tetracc.c -o $@

#-----------------------------------------------------------------------
# Utility Functions
#-----------------------------------------------------------------------
binarize.o: $(HDRS)
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
        zip -rq corr.zip    corr/{src,doc} cpuinfo/{src,doc}; \
        tar cfz corr.tar.gz corr/{src,doc} cpuinfo/{src,doc}

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
