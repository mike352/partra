staticlib: objects

objects: ising.o potts.o reductions.o genfuncs.o
	ar rcs libpartra.a ising.o potts.o reductions.o genfuncs.o

ising.o: ising.c
	gcc -c -Ofast ising.c

potts.o: potts.c
	gcc -c -Ofast potts.c

reductions.o: reductions.c
	gcc -c -Ofast reductions.c

genfuncs.o: genfuncs.c
	gcc -c -Ofast genfuncs.c

clean: 
	rm *.o staticlib