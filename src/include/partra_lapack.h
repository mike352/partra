#ifndef PARTRA_LAPACK_H
#define PARTRA_LAPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>
#include "lapacke.h"

unsigned char eigs_lapack(double****, unsigned long long*); 

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
#endif

#endif