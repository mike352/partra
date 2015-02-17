# Compiler to use:
CC=gcc

#Flags to use when compiling:
CFLAGS1=-c -Wall
#-O3 is safe. -Ofast should be tested for numerics, but safe for integers
CFLAGS2=-Ofast

staticlib: objects

objects: ising.o potts.o reductions.o genfuncs.o
	ar rcs libpartra.a ising.o potts.o reductions.o genfuncs.o

ising.o: ising.c
	$(CC) $(CFLAGS1) $(CFLAGS2) ising.c

potts.o: potts.c
	$(CC) $(CFLAGS1) $(CFLAGS2) potts.c

reductions.o: reductions.c
	$(CC) $(CFLAGS1) $(CFLAGS2) reductions.c

genfuncs.o: genfuncs.c
	$(CC) $(CFLAGS1) $(CFLAGS2) genfuncs.c

clean: 
	rm *.o staticlib