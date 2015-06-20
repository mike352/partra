#include "partra_genfuncs.h"
#include "partra_ising.h"

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

/*****************************************************/
/**********Full Ising square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_sq_f_f_%d",N);

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
		(*matrix)[n][m][1][0]=bit_sum(n^m)+uh;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char i_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_sq_c_f_%d",N);

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
		(*matrix)[n][m][1][0]=bit_sum(n^m)+uh;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_sq_f_f_%d",N);

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
		(*matrix)[n][m][1][0]=bit_sum(n^m)+uh;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_sq_c_f_%d",N);

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
		(*matrix)[n][m][1][0]=bit_sum(n^m)+uh;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*****************************************************/
/*******************SYMMETRIC MATRICES****************/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_f_s(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_sq_f_f_s_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,uh2;
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
		uh2 = bit_sum((~1ULL&m)^(~1ULL&circ_single_lshift(m,N)));
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=2*bit_sum(n^m)+uh+uh2;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char i_sq_c_f_s(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
sprintf(filename,"i_sq_c_f_s_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,uh2;
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
		uh2 = bit_sum(m^circ_single_lshift(m,N));
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=2*bit_sum(n^m)+uh+uh2;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_sq_f_f_s(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_sq_f_f_s_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,uh2,xh,xh2;
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
		uh2 = bit_sum((~1ULL&m)^(~1ULL&circ_single_lshift(m,N)));
		xh2 = N-(bit_sum(m));
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=2*bit_sum(n^m)+uh+uh2;
		(*matrix)[n][m][1][1]=xh+xh2;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*******************************/
unsigned char if_sq_c_f_s(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
sprintf(filename,"if_sq_c_f_s_%d",N);

unsigned long long n;
unsigned long long m;
unsigned char uh,uh2,xh,xh2;
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
		uh2 = bit_sum(m^circ_single_lshift(m,N));
		xh2 = N-(bit_sum(m));
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=2*bit_sum(n^m)+uh+uh2;
		(*matrix)[n][m][1][1]=xh+xh2;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}
