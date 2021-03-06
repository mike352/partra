#ifndef PARTRA_LAPACK_H
#define PARTRA_LAPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>

unsigned char eigs_lapack_file(double****, unsigned long long*,char*); 

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
typedef double complex**** partra_matrix_dc;
#endif

#endif