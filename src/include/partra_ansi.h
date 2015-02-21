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

//Partra header files
#include "partra_genfuncs.h"
#include "partra_reductions.h"
#include "partra_ising.h"
#include "partra_ising_ldc.h"
#include "partra_potts.h"

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** Matrix;
typedef unsigned long long**** Matrix_ll;
typedef long double complex**** Matrix_ldc;
typedef unsigned char*** Row;
typedef unsigned long long*** Row_ll;
typedef long double complex*** Row_ldc;
#endif

#endif