#ifndef PARTRA_GENFUNCS_D_H
#define PARTRA_GENFUNCS_D_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

unsigned char matrix_alloc_d(double*****,const unsigned long long*, const unsigned char);
void matrix_free_d(double****, const unsigned long long*);

unsigned char matrix_fprintf_d(double****, unsigned long long*, char*,unsigned char);

unsigned char matrix_sub_d(double*****, unsigned long long*, unsigned char****, unsigned long long*, char*, ...);
unsigned char matrix_sub_d_d(double*****, unsigned long long*, double****, unsigned long long*, char*, double);


//Double Complex Functions
/*
unsigned char matrix_alloc_dc(double complex*****,const unsigned long long*, const unsigned char);
void matrix_free_dc(double complex****, const unsigned long long*);

unsigned char matrix_fprintf_dc(double complex****, unsigned long long*, char*,unsigned char);

unsigned char matrix_sub_dc(double complex*****, unsigned long long*, unsigned char****, unsigned long long*, char*, ...);
unsigned char matrix_sub_dc_dc(double complex*****, unsigned long long*, double complex****, unsigned long long*, double complex);
unsigned char matrix_sub_dc_d(double complex*****, unsigned long long*, double****, unsigned long long*, double complex);
*/
#ifdef __cplusplus
}
#endif

#endif