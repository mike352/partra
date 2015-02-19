# Compiler to use:
CC=gcc

#Flags to use when compiling:
CFLAGS1=-Wall
CFLAGS2=-c
#Optimization flags: -O3 is safe. -Ofast should be tested for numerics, but safe for integers
COFLAG=-Ofast

OBJFILES= $(addsuffix .o,$(basename $(wildcard src/ising/*.c))) \
		  $(addsuffix .o,$(basename $(wildcard src/potts/*.c)))  \
		  $(addsuffix .o,$(basename $(wildcard src/reductions/*.c))) \
		  $(addsuffix .o,$(basename $(wildcard src/genfuncs/*.c)))

#Header files
HEADFILES := $(wildcard src/include/*.h)

#Header directory
HEADDIR=src/include

#Static library location
LIBDIR=src/lib

#Path to place executibles, headers, and static libraries
PREFIX=/usr


compile: staticlib standalone

standalone: staticlib partra examples

partra: src/partra.c
	$(CC) $(CFLAGS1) $(COFLAG) src/partra.c -lpartra -I$(HEADDIR) -L. -o partra

examples: bitarch example

bitarch: examples/bitarchitecture.c
	$(CC) $(CFLAGS1) $(COFLAG) examples/bitarchitecture.c -lpartra  -I$(HEADDIR) -L. -o partra_bitarch

example: examples/example.c
	$(CC) $(CFLAGS1) $(COFLAG) examples/example.c -lpartra  -I$(HEADDIR) -L. -o partra_example
	
staticlib: objects

objects: $(OBJFILES)
	ar rcs libpartra.a $^

src/ising/%.o: src/ising/%.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) $< -o $@ -I$(HEADDIR)

src/potts/%.o: src/potts/%.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) $< -o $@ -I$(HEADDIR)

src/reductions/%.o: src/reductions/%.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) $< -o $@ -I$(HEADDIR)

src/genfuncs/%.o: src/genfuncs/%.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(COFLAG) $< -o $@ -I$(HEADDIR)

install:
	cp libpartra.a $(LIBDIR)/
	mv libpartra.a $(PREFIX)/lib
	cp $(HEADFILES) $(PREFIX)/include
	mv partra partra_bitarch partra_example $(PREFIX)/bin

all: compile install
	
clean: 
	rm $(OBJFILES)