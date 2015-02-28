#ifndef PARTRA_H
#define PARTRA_H

#include "partra_ansi.h"
#include "partra_gmp.h"
#include "partra_lapack.h"

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
#endif

#endif