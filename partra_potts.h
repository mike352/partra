#ifndef PARTRA_POTTS_H
#define PARTRA_POTTS_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

unsigned char p_sq_f_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_tri_f_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p2_sq_f_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

unsigned char p2_tri_f_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

#ifdef __cplusplus
}
#endif

#endif