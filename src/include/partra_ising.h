#ifndef PARTRA_ISING_H
#define PARTRA_ISING_H

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
//typedef double complex**** partra_matrix_dc;
#endif

#ifdef __cplusplus
extern "C" {
#endif

//Full and reduced transfer matrices, square lattice, for creating files
unsigned char i_sq_f_f_file(const unsigned char, const char*); //Ising full transfer matrix, free row b.c.
unsigned char i_sq_c_f_file(const unsigned char, const char*); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_sq_f_f_file(const unsigned char, const char*); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_sq_c_f_file(const unsigned char, const char*); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_sq_f_r_file(const unsigned char, const char*); //Ising reduced transfer matrix, free row b.c.
unsigned char i_sq_c_r_file(const unsigned char, const char*); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_sq_f_r_file(const unsigned char, const char*); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_sq_c_r_file(const unsigned char, const char*); //Ising in a field reduced transfer matrix, cylindrical row b.c.

//Full transfer matrices, square lattice, for creating files, symmetric matrices
unsigned char i_sq_f_f_s_file(const unsigned char, const char*); //Ising full transfer matrix, free row b.c., symmetric matrix
unsigned char i_sq_c_f_s_file(const unsigned char, const char*); //Ising full transfer matrix, cylindrical row b.c., symmetric matrix
unsigned char if_sq_f_f_s_file(const unsigned char, const char*); //Ising in a field full transfer matrix, free row b.c., symmetric matrix
unsigned char if_sq_c_f_s_file(const unsigned char, const char*); //Ising in a field full transfer matrix, cylindrical row b.c., symmetric matrix

//Full and reduced transfer matrices, triangular lattice, for creating files
unsigned char i_tri_f_f_file(const unsigned char, const char*); //Ising full transfer matrix, free row b.c.
unsigned char i_tri_c_f_file(const unsigned char, const char*); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_tri_f_f_file(const unsigned char, const char*); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_tri_c_f_file(const unsigned char, const char*); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_tri_f_r_file(const unsigned char, const char*); //Ising reduced transfer matrix, free row b.c.
unsigned char i_tri_c_r_file(const unsigned char, const char*); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_tri_f_r_file(const unsigned char, const char*); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_tri_c_r_file(const unsigned char, const char*); //Ising in a field reduced transfer matrix, cylindrical row b.c.

//Full and reduced transfer matrices, square lattice
unsigned char i_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising full transfer matrix, free row b.c.
unsigned char i_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_sq_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_sq_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising reduced transfer matrix, free row b.c.
unsigned char i_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_sq_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_sq_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field reduced transfer matrix, cylindrical row b.c.

//Full transfer matrices, triangular lattice
unsigned char i_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising full transfer matrix, free row b.c.
unsigned char i_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_tri_f_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_tri_c_f(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field full transfer matrix, cylindrical row b.c.
/*
unsigned char i_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising reduced transfer matrix, free row b.c.
unsigned char i_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_tri_f_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_tri_c_r(unsigned char*****, unsigned long long*, char*, const unsigned char); //Ising in a field reduced transfer matrix, cylindrical row b.c.
*/
#ifdef __cplusplus
}
#endif

#endif