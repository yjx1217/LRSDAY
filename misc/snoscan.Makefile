# Makefile for snoscan
# TML, 12/8/96
#
# Makes Sean Eddy's sequence handling library (squid),
# then compiles snoscan program using squid library.
#

VERSION = "0.9.1"
VDATE   = "26 Apr 2017"

## where you want things installed
# BINDIR  = $(HOME)/bin
# MANDIR  = $(HOME)/man

## where are the squid library and include files?
SQUIDDIR = ./squid-1.5.11
## and the squid source code?
SQUIDSRCDIR = ./squid-1.5.11

## your compiler  -- must be same as Makefile in "squid" directory!
CC = cc         # use cc if you don't have gcc
#CC = gcc	# use GNU gcc, if you have it...

## any special compiler flags you want
##  -- must be same as Makefile in "squid" directory!
CFLAGS = -O2                            # most machines, production build
#CFLAGS = -g                            # gcc, development/maintenance



## machine specific definitions
# -DNO_STRDUP    this machine has no strdup() function
# -DNOSTR        this machine has no strstr() function
# -DNORANDOM     this machine has no random() function (Solaris)


# Most machines need nothing (includes Sun SPARC)
MDEFS = 
# Uncomment this one for debugging version, Silicon Graphics
#MDEFS = -DNEED_GETOPTH -DMEMDEBUG -DDEBUG


#######
## should not need to modify below this line
#######
SHELL  = /bin/sh
LIBS   = -lsquid -lm

VFLAGS  = -DVERSION=$(VERSION) -DVDATE=$(VDATE)

READMES = README 

PROGS = snoscan

HDRS =  snoscan.h matrices.h

MAIN =  snoscan_main.c

SQUIDSRC = alignio.c sqerror.c sqio.c iupac.c msf.c revcomp.c\
	selex.c sre_ctype.c sre_string.c stack.c types.c

SQUIDHDRS = squid.h sqfuncs.h

SRC =   search.c

OBJ =   search.o

DISTFILES = $(MAIN) $(SRC) $(HDRS) $(MANSRC) $(READMES) $(DOCS)


all: 	$(PROGS)

snoscan: $(OBJ) snoscan_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o snoscan -L$(SQUIDDIR) snoscan_main.o $(OBJ) $(LIBS)
	chmod +x sort-snos scan-yeast

install: $(PROGS)
	cp $(PROGS) sort-snos $(BINDIR)/

clean:
	rm -f *.o *~ core $(PROGS) TAGS *.raw *.sort-all *.sort-bysite *.sortRH

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<		

.m.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<





