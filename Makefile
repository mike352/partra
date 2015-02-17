# Compiler to use:
CC=gcc

#Flags to use when compiling:
CFLAGS1=-Wall
CFLAGS2=-c 
#Optimization flags: -O3 is safe. -Ofast should be tested for numerics, but safe for integers
COFLAG=-Ofast

#Linking flags of the newly compiled library
CLFLAGS=-lpartra -L.

all: staticlib standalone

standalone: staticlib partra

partra: partra.c
	$(CC) $(CFLAGS1) $(COFLAG) partra.c $(CLFLAGS)

staticlib: objects

objects: ising.o potts.o reductions.o genfuncs.o
	ar rcs libpartra.a ising.o potts.o reductions.o genfuncs.o

ising.o: ising.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) ising.c

potts.o: potts.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) potts.c

reductions.o: reductions.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) reductions.c

genfuncs.o: genfuncs.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) genfuncs.c

clean: 
	rm *.o staticlib