#-----------------------------------------------------------------------
# File    : makefile
# Contents: build pcc and tetracc programs (on Unix systems)
# Author  : Christian Borgelt
#-----------------------------------------------------------------------
SHELL       = /bin/bash
CPUINFODIR  = ../../cpuinfo/cpuinfo/src
BINARIZEDIR = ../../binarize/binarize/src

CC       = gcc -std=c99 -march=native
# CC       = gcc -std=c99 -mavx
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3
# CFLAGS   = $(CFBASE) -g
INCS     = -I$(CPUINFODIR) -I$(BINARIZEDIR)

LD       = gcc
LDFLAGS  =
LIBS     = -lm -lrt -lpthread

HDRS     = $(BINARIZEDIR)/binarize.h $(CPUINFODIR)/cpuinfo.h
OBJS     = $(BINARIZEDIR)/binarize.o $(CPUINFODIR)/cpuinfo.o
PRGS     = pcc tetracc

#-----------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------
all:        $(PRGS)

pcc:        $(OBJS) pcc.o makefile
	$(LD) $(LDFLAGS) $(OBJS) pcc.o     $(LIBS) -o $@

tetracc:    $(OBJS) tetracc.o makefile
	$(LD) $(LDFLAGS) $(OBJS) tetracc.o $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
pcc.o:      pcc.h $(HDRS)
pcc.o:      pcc.c makefile
	$(CC) $(CFLAGS) $(INCS) -DPCC_MAIN -c pcc.c     -o $@

tetracc.o:  tetracc.h $(HDRS)
tetracc.o:  tetracc.c tetracc.h makefile
	$(CC) $(CFLAGS) $(INCS) -DTCC_MAIN -c tetracc.c -o $@

#-----------------------------------------------------------------------
# Utility Functions
#-----------------------------------------------------------------------
binarize.o: $(HDRS)
binarize.o: binarize.c makefile
	$(CC) $(CFLAGS) $(INCS) -c binarize.c -o $@

cpuinfo.o:  $(HDRS)
cpuinfo.o:  cpuinfo.c makefile
	$(CC) $(CFLAGS) $(INCS) -c cpuinfo.c -o $@

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f pcc.zip pcc.tar.gz; \
        zip -rq pcc.zip    pcc/{src,ex,doc}; \
        tar cfz pcc.tar.gz pcc/{src,ex,doc}

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) $(HOME)/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.o *~ *.flc core $(PRGS)

clean:
	$(MAKE) localclean
