ARCH ?= $(shell uname -m)

CC = gcc
GCC = g++
CDEBUG =#-ggdb -g -pg
OPTM =-O3
WARN = -W
KLIB = ../lib/kent/src/lib/$(ARCH)
KINC = ../lib/kent/src/inc
CFLAGS = $(WARN) $(OPTM) -I. -I$(KINC)
CLIB = $(KLIB)/jkweb.a -lm

BIN = $(HOME)/bin/$(ARCH)
RM = rm -rf

ALLSRC = inferAdjProb raca

all: $(ALLSRC)

%: %.c
	$(CC) $(CDEBUG) $(CFLAGS) \
   $+ \
   $(CLIB) \
  -o $@

%: %.cpp
	$(GCC) \
   $+ \
  -o $@

.PHONY: tags
tags:
	ctags *.[hc] lib/*.[hc]	inc/*.h

.PHONY: clean
clean:
	$(RM) $(ALLSRC) *.o *.dSYM
	
