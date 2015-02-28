#ifndef PARTRA_GENFUNCS_MPF_H
#define PARTRA_GENFUNCS_MPF_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <complex.h>
#include <stdarg.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
	
unsigned char matrix_sub_d_mpf(double*****, unsigned long long*, unsigned char****, unsigned long long*, char*, ...);

unsigned char matrix_sub_d_d_mpf(double*****, unsigned long long*, double****, unsigned long long*, double);

#ifdef __cplusplus
}
#endif

#endif