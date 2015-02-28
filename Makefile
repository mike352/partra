# Compiler to use:
CC=gcc

#Flags to use when compiling:
CFLAGS=-Wall -lm 
CFLAGSo=-c
LIB1=-lgmp
LIB2=-lblas
LIB3=-llapack
LIB4=-larpack
LIB5=-lmpc
#Optimization flags: -O3 is safe. -Ofast should be tested for numerics, but safe for integers
COFLAG=-Ofast

OBJFILES= $(addsuffix .o,$(basename $(wildcard src/ising/*.c))) \
	  $(addsuffix .o,$(basename $(wildcard src/potts/*.c)))  \
	  $(addsuffix .o,$(basename $(wildcard src/reductions/*.c))) \
	  $(addsuffix .o,$(basename $(wildcard src/genfuncs/*.c))) \
	  $(addsuffix .o,$(basename $(wildcard src/eigenvalues/*.c)))

EXPOUT= $(addprefix partra_, $(notdir $(basename $(wildcard examples/*.c))))
		  
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

partra: src/partra.c libpartra.a
	$(CC) src/partra.c -lpartra -I$(HEADDIR) -L. -o partra $(CFLAGS) $(COFLAG)

examples: $(EXPOUT)

partra_%: examples/%.c
	$(CC) $< -lpartra  -I$(HEADDIR) -L. $(addprefix -o , $@) $(CFLAGS) $(COFLAG)
	
staticlib: libpartra.a

libpartra.a: $(OBJFILES)
	ar rcs libpartra.a $^

src/ising/%.o: src/ising/%.c
	$(CC) $< -o $@ -I$(HEADDIR) $(CFLAGSo) $(CFLAGS) $(COFLAG)

src/potts/%.o: src/potts/%.c
	$(CC) $< -o $@ -I$(HEADDIR) $(CFLAGSo) $(CFLAGS) $(COFLAG)

src/reductions/%.o: src/reductions/%.c
	$(CC) $< -o $@ -I$(HEADDIR) $(CFLAGSo) $(CFLAGS) $(COFLAG)

src/genfuncs/%.o: src/genfuncs/%.c
	$(CC) $< -o $@ -I$(HEADDIR) $(CFLAGSo) $(PACK1) $(CFLAGS) $(COFLAG)

src/eigenvalues/%.o: src/eigenvalues/%.c
	$(CC) $< -o $@ -I$(HEADDIR) $(CFLAGSo) $(PACK2) $(PACK3) $(CFLAGS) $(COFLAG)

install:
	cp -a libpartra.a $(LIBDIR)/
	mv libpartra.a $(PREFIX)/lib
	cp -a $(HEADFILES) $(PREFIX)/include
	mv partra $(PREFIX)/bin
	mv $(EXPOUT) $(PREFIX)/bin

all: compile install
	
clean: 
	rm $(OBJFILES)