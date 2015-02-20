#ifndef PARTRA_GENFUNCS_H
#define PARTRA_GENFUNCS_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <complex.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif
	
/*This function performs a circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_single_lshift(unsigned long long, const unsigned char);

/*This function performs a binned circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_bin_lshift(unsigned long long, const unsigned char, const unsigned char);

/*This function reflects all the bits for a number of a particular width*/
unsigned long long bit_reflection(unsigned long long, const unsigned char);

/*This function reflects all the bits for a number of a particular width and bin size bin*/
unsigned long long bit_reflection_bin(unsigned long long, const unsigned char, const unsigned char);

/*This function adds up all of the set bits in a number x*/
unsigned char bit_sum(unsigned long long);


/*A not very efficient function for taking the power of a matrix*/
unsigned char matrix_pow_ll(unsigned long long*****, unsigned long long*****, const unsigned long long*,const unsigned long long);

/*Matrix and row allocation routines*/
unsigned char matrix_alloc(unsigned char*****,const unsigned long long*,const unsigned char);
void matrix_free(unsigned char****,const unsigned long long*);
unsigned char row_alloc(unsigned char****,const unsigned long long*,const unsigned char);
void row_free(unsigned char***,const unsigned long long*);
unsigned char matrix_setadd(unsigned char*****, const unsigned long long*, const unsigned long long, const unsigned long long, const unsigned char*);
unsigned char row_setadd(unsigned char****, const unsigned long long*, const unsigned long long, const unsigned char*);

unsigned char matrix_alloc_ll(unsigned long long*****,const unsigned long long*,const unsigned char);
void matrix_free_ll(unsigned long long****,const unsigned long long*);
unsigned char row_alloc_ll(unsigned long long****,const unsigned long long*,const unsigned char);
void row_free_ll(unsigned long long***,const unsigned long long*);
unsigned char matrix_setadd_ll(unsigned long long*****, const unsigned long long*, const unsigned long long, const unsigned long long, const unsigned long long*);
unsigned char row_setadd_ll(unsigned long long****, const unsigned long long*, const unsigned long long, const unsigned long long*);

/*Matrix functions*/
unsigned char matrix_fprintf(const unsigned char*****, const unsigned long long*, const char*, const char*);

#ifdef __cplusplus
}
#endif

#endif