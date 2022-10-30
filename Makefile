SHELL=/bin/sh
CC=gcc 
CFLAGS=-Wall -pedantic -ansi -O3 -D_GNU_SOURCE
AR=ar

TARGETS=optsn

.PHONY: all clean

all: $(TARGETS)

lmfit.o: lmfit.c lmfit.h	
	$(CC) $(CFLAGS) -c lmfit.c

downhill.o: downhill.c downhill.h
	$(CC) $(CFLAGS) -c downhill.c 

optsn: optsn.c lmfit.o downhill.o
	$(CC) $(CFLAGS) -o optsn optsn.c lmfit.o downhill.o -lm

clean:
	rm -f $(TARGETS) *.o
