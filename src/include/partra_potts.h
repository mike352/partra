#ifndef PARTRA_POTTS_H
#define PARTRA_POTTS_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <complex.h>

#ifndef PARTRA_TYPEDEF
#define PARTRA_TYPEDEF
typedef unsigned char**** partra_matrix;
typedef unsigned char*** partra_row;

typedef double**** partra_matrix_d;
typedef double complex**** partra_matrix_dc;
#endif

#ifdef __cplusplus
extern "C" {
#endif

unsigned char p_sq_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_sq_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_tri_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_tri_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p2_sq_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_sq_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

unsigned char p2_tri_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_f_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_tri_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2


unsigned char p_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c.
/*
unsigned char p_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c.
*/
unsigned char p_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c.
/*
unsigned char p_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c.
*/
unsigned char p2_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
/*
unsigned char p2_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2
*/
unsigned char p2_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
/*
unsigned char p2_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2
*/



//Symmetric matrices
unsigned char p_sq_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_sq_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_tri_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_tri_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p2_sq_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_sq_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

unsigned char p2_tri_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_f_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_tri_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r_s_file(const unsigned char, const unsigned long long, const char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2


unsigned char p_sq_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c.
/*
unsigned char p_sq_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c.
*/
unsigned char p_tri_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c.
/*
unsigned char p_tri_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c.
*/
unsigned char p2_sq_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
/*
unsigned char p2_sq_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2
*/
unsigned char p2_tri_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_f_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
/*
unsigned char p2_tri_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r_s(unsigned char*****, unsigned long long*, char*, const unsigned char, const unsigned long long); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2
*/




#ifdef __cplusplus
}
#endif

#endif