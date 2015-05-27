#ifndef PARTRA_ANSI_H
#define PARTRA_ANSI_H

#include "partra_genfuncs.h"
#include "partra_genfuncs_d.h"
#include "partra_reductions.h"
#include "partra_ising.h"
//#include "partra_ising_ldc.h"
#include "partra_potts.h"

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
typedef double complex**** partra_matrix_dc;
#endif

#endif