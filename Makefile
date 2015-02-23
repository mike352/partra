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
	$(CC) $(CFLAGS1) $(COFLAG) src/partra.c -lpartra -I$(HEADDIR) -L. -o partra

examples: $(EXPOUT)

partra_%: examples/%.c
	$(CC) $(CFLAGS1) $(COFLAG) $< -lpartra  -I$(HEADDIR) -L. $(addprefix -o , $@)
	
staticlib: libpartra.a

libpartra.a: $(OBJFILES)
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
	cp -a libpartra.a $(LIBDIR)/
	mv libpartra.a $(PREFIX)/lib
	cp -a $(HEADFILES) $(PREFIX)/include
	mv partra $(PREFIX)/bin
	mv $(EXPOUT) $(PREFIX)/bin

all: compile install
	
clean: 
	rm $(OBJFILES)