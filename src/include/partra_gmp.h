#ifndef PARTRA_ANSI_H
#define PARTRA_ANSI_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <complex.h>
#include <time.h>
#include "gmp.h"

//Partra header files
#include "partra_genfuncs_mpf.h"

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
#endif

#endif