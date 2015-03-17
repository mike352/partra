#include "partra_genfuncs.h"
#include "partra_potts.h"


/*****************************************************/
/*********Full Potts square transfer matrices*********/
/*****************************************************/

/*******************************/
unsigned char p2_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=2ULL;
sprintf(filename,"p_sq_f_f_%llu_%d",Q,N);

unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	utest = n ^ circ_bin_lshift(n,N,bin);
	uh=0ULL;
	for (p=1;p<N;p++)
	{
		uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}

	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{	
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=uh+uinter;
		(*matrix)[n][m][1][1]=1;
	}
}
	
return 0;
}


/*******************************/
unsigned char p2_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=2ULL;
sprintf(filename,"p_sq_c_f_%llu_%d",Q,N);

unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	utest = n ^ circ_bin_lshift(n,N,bin);
	uh=0ULL;
	for (p=0;p<N;p++)
	{
		uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=uh+uinter;
		(*matrix)[n][m][1][1]=1;
	}
}

return 0;
}


/*******************************/
unsigned char pf2_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=3ULL;
sprintf(filename,"pf_sq_f_f_%llu_%d",Q,N);

unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter,xh;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	utest = n ^ circ_bin_lshift(n,N,bin);
	uh=0ULL;
	for (p=1;p<N;p++)
	{
		uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	xh = 0ULL;
	for (p=0;p<N;p++)
	{
		xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=uh+uinter;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*******************************/
unsigned char pf2_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=3ULL;
sprintf(filename,"pf_sq_c_f_%llu_%d",Q,N);

unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter,xh;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	utest = n ^ circ_bin_lshift(n,N,bin);
	uh=0ULL;
	for (p=0;p<N;p++)
	{
		uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	xh = 0ULL;
	for (p=0;p<N;p++)
	{
		xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		(*matrix)[n][m][0][0]=1;
		(*matrix)[n][m][1][0]=uh+uinter;
		(*matrix)[n][m][1][1]=xh;
		(*matrix)[n][m][1][2]=1;
	}
}

return 0;
}


/*******************************/
unsigned char p_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=2ULL;
sprintf(filename,"p_sq_f_f_%llu_%d",Q,N);

const unsigned char csize=CHAR_BIT;
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	matrix_free(*matrix,msize);
	return 2;
}	

//Find valid Potts numbers
for (p=0ULL;p<(1ULL<<bin)-Q;p++)
{	
	for (n=0ULL;n<(1ULL<<(bin*N));n++)
	{
		sum = 0ULL;
		for (m=0ULL;m<N;m++)
		{
			sum = sum + (((n&(((1ULL<<bin)-p-1)<<bin*m))>>bin*m)==((1ULL<<bin)-p-1));
		}
		if (sum!=0)
		{
			bitfrac=lldiv(n,csize); 
			pnums[bitfrac.quot]=pnums[bitfrac.quot] | (1<<bitfrac.rem);
		}
	}
}
//Convention: pnums is a bitarray whose 0 bits correspond to the valid Potts numbers

//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	bitfrac=lldiv(n,csize);
	if (((pnums[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				(*matrix)[n][m][0][0]=1;
				(*matrix)[n][m][1][0]=uh+uinter;
				(*matrix)[n][m][1][1]=1;
			}
		}
	}
}
	

free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=2ULL;
sprintf(filename,"p_sq_c_f_%llu_%d",Q,N);

const unsigned char csize=CHAR_BIT;
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	matrix_free(*matrix,msize);
	return 2;
}	
	
//Find valid Potts numbers
for (p=0ULL;p<(1ULL<<bin)-Q;p++)
{	
	for (n=0ULL;n<(1ULL<<(bin*N));n++)
	{
		sum = 0ULL;
		for (m=0ULL;m<N;m++)
		{
			sum = sum + (((n&(((1ULL<<bin)-p-1)<<bin*m))>>bin*m)==((1ULL<<bin)-p-1));
		}
		if (sum!=0)
		{
			bitfrac=lldiv(n,csize); 
			pnums[bitfrac.quot]=pnums[bitfrac.quot] | (1<<bitfrac.rem);
		}
	}
}	
//Convention: pnums is a bitarray whose 0 bits correspond to the valid Potts numbers

//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	bitfrac=lldiv(n,csize);
	if (((pnums[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				(*matrix)[n][m][0][0]=1;
				(*matrix)[n][m][1][0]=uh+uinter;
				(*matrix)[n][m][1][1]=1;
			}
		}
	}
}

free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_f_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=3ULL;
sprintf(filename,"pf_sq_f_f_%llu_%d",Q,N);

const unsigned char csize=CHAR_BIT;
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter,xh;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	matrix_free(*matrix,msize);
	return 2;
}	

//Find valid Potts numbers
for (p=0ULL;p<(1ULL<<bin)-Q;p++)
{	
	for (n=0ULL;n<(1ULL<<(bin*N));n++)
	{
		sum = 0ULL;
		for (m=0ULL;m<N;m++)
		{
			sum = sum + (((n&(((1ULL<<bin)-p-1)<<bin*m))>>bin*m)==((1ULL<<bin)-p-1));
		}
		if (sum!=0)
		{
			bitfrac=lldiv(n,csize); 
			pnums[bitfrac.quot]=pnums[bitfrac.quot] | (1<<bitfrac.rem);
		}
	}
}	
//Convention: pnums is a bitarray whose 0 bits correspond to the valid Potts numbers

//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	bitfrac=lldiv(n,csize);
	if (((pnums[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh = 0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				(*matrix)[n][m][0][0]=1;
				(*matrix)[n][m][1][0]=uh+uinter;
				(*matrix)[n][m][1][1]=xh;
				(*matrix)[n][m][1][2]=1;
			}
		}
	}
}

free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_c_f(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N, const unsigned long long Q)
{
msize[1]=3ULL;
sprintf(filename,"pf_sq_c_f_%llu_%d",Q,N);

const unsigned char csize=CHAR_BIT;
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter,xh;
unsigned char flag;

while((1ULL<<bin)<Q)
{
	bin++;
}

msize[0] = 1ULL<<N;
flag = matrix_alloc(matrix,msize,1);
if (flag!=0)
{
	return flag;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	matrix_free(*matrix,msize);
	return 2;
}	

//Find valid Potts numbers
for (p=0ULL;p<(1ULL<<bin)-Q;p++)
{	
	for (n=0ULL;n<(1ULL<<(bin*N));n++)
	{
		sum = 0ULL;
		for (m=0ULL;m<N;m++)
		{
			sum = sum + (((n&(((1ULL<<bin)-p-1)<<bin*m))>>bin*m)==((1ULL<<bin)-p-1));
		}
		if (sum!=0)
		{
			bitfrac=lldiv(n,csize); 
			pnums[bitfrac.quot]=pnums[bitfrac.quot] | (1<<bitfrac.rem);
		}
	}
}	
//Convention: pnums is a bitarray whose 0 bits correspond to the valid Potts numbers

//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	bitfrac=lldiv(n,csize);
	if (((pnums[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh = 0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				(*matrix)[n][m][0][0]=1;
				(*matrix)[n][m][1][0]=uh+uinter;
				(*matrix)[n][m][1][1]=xh;
				(*matrix)[n][m][1][2]=1;
			}
		}
	}
}

free((void*)pnums);
return 0;
}

