#ifndef PARTRA_REDUCTIONS_H
#define PARTRA_REDUCTIONS_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

//Partra header files
#include "partra_genfuncs.h"

extern "C"
{
	
unsigned char red_simple_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N
unsigned char red_simple_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N
unsigned char red_simple_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N*bin
unsigned char red_simple_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N*bin

unsigned char red_gen_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array

}

#endif