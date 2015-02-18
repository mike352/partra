# Compiler to use:
CC=gcc

#Flags to use when compiling:
CFLAGS1=-Wall
CFLAGS2=-c
#Optimization flags: -O3 is safe. -Ofast should be tested for numerics, but safe for integers
COFLAG=-Ofast

#Linking flags of the newly compiled library
CLFLAGS=-lpartra -L.

#Change this path if placing executibles, headers, and static libraries in a different location
PREFIX=/usr


compile: staticlib standalone

standalone: staticlib partra example bitarch

partra: partra.c
	$(CC) $(CFLAGS1) $(COFLAG) partra.c $(CLFLAGS) -o partra

example: example.c
	$(CC) $(CFLAGS1) $(COFLAG) example.c $(CLFLAGS) -o partra_example

bitarch: bitarchitecture.c
	$(CC) $(CFLAGS1) $(COFLAG) bitarchitecture.c $(CLFLAGS) -o partra_bitarch
	
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

install:
	mv libpartra.a $(PREFIX)/lib
	cp *.h $(PREFIX)/include
	mv partra partra_bitarch partra_example $(PREFIX)/bin

all: compile install
	
clean: 
	rm *.o