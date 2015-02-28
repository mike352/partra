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

#ifdef __cplusplus
extern "C" {
#endif

unsigned char matrix_alloc_d(double*****,const unsigned long long*, const unsigned char);
void matrix_free_d(double****, const unsigned long long*);

/*Matrix functions*/
unsigned char matrix_fprintf_d(double****, unsigned long long*, char*,unsigned char);
unsigned char matrix_sub_d(double*****, unsigned long long*, unsigned char****, unsigned long long*, char*, ...);
unsigned char matrix_sub_d_d(double*****, unsigned long long*, double****, unsigned long long*, double);

#ifdef __cplusplus
}
#endif

#endif