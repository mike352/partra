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
#include "gmp.h"

#ifdef __cplusplus
extern "C" {
#endif
	
unsigned char matrix_sub_mpf_d(double*****, unsigned long long*, unsigned char****, unsigned long long*, char*, ...);

unsigned char matrix_sub_mpf_d_d(double*****, unsigned long long*, double****, unsigned long long*, double);

#ifdef __cplusplus
}
#endif

#endif