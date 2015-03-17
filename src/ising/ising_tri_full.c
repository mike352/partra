#include "partra_genfuncs.h"
#include "partra_ising.h"


/*****************************************************/
/*******Full Ising triangular transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char i_tri_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_tri_f_f_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh;
unsigned char flag;

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	for(m=0; m<(1ULL<<N); m++)
	{
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char i_tri_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_tri_c_f_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh;
unsigned char flag;

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum(n^circ_single_lshift(n,N));
	for(m=0; m<(1ULL<<N); m++)
	{
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_tri_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_tri_f_f_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,xh;
unsigned char flag;

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	xh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_tri_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_tri_c_f_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,xh;
unsigned char flag;

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum(n^circ_single_lshift(n,N));
	xh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}
