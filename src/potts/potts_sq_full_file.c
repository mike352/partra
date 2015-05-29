#include "partra_genfuncs.h"
#include "partra_potts.h"


/*****************************************************/
/*********Full Potts square transfer matrices*********/
/*****************************************************/

/*******************************/
unsigned char p2_sq_f_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
	
		fprintf(fid,"%hhu\n",uh+uinter);
	}
}
	
printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p2_sq_c_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		
		fprintf(fid,"%hhu\n",uh+uinter);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_f_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		
		fprintf(fid,"%hhu %hhu\n",uh+uinter,xh);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_c_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		
		fprintf(fid,"%hhu %hhu\n",uh+uinter,xh);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p_sq_f_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
			
				fprintf(fid,"%hhu\n",uh+uinter);
			}
		}
	}
}
	

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_sq_c_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
				
				fprintf(fid,"%hhu\n",uh+uinter);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_f_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
				
				fprintf(fid,"%hhu %hhu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_c_f_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
				
				fprintf(fid,"%hhu %hhu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}





/*****************************************************/
/******************SYMMETRIC MATRICES*****************/
/*****************************************************/

/*******************************/
unsigned char p2_sq_f_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uh2,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		utest = m ^ circ_bin_lshift(m,N,bin);
		uh2=0ULL;
		for (p=1;p<N;p++)
		{
			uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	
		fprintf(fid,"%hhu\n",uh+uh2+2*uinter);
	}
}
	
printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p2_sq_c_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uh2,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		utest = m ^ circ_bin_lshift(m,N,bin);
		uh2=0ULL;
		for (p=0;p<N;p++)
		{
			uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	  
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%hhu\n",uh+uh2+2*uinter);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_f_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uh2,uinter,xh,xh2;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		utest = m ^ circ_bin_lshift(m,N,bin);
		uh2=0ULL;
		for (p=1;p<N;p++)
		{
			uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh2 = 0ULL;
		for (p=0;p<N;p++)
		{
			xh2 = xh2 + (((m & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	  
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%hhu %hhu\n",uh+uh2+2*uinter,xh+xh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_c_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest;
unsigned char uh,uh2,uinter,xh,xh2;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
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
		utest = m ^ circ_bin_lshift(m,N,bin);
		uh2=0ULL;
		for (p=0;p<N;p++)
		{
			uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh2 = 0ULL;
		for (p=0;p<N;p++)
		{
			xh2 = xh2 + (((m & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	  
		utest = n ^ m;
		uinter=0ULL;
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%hhu %hhu\n",uh+uh2+2*uinter,xh+xh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p_sq_f_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uh2,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
			utest = m ^ circ_bin_lshift(m,N,bin);
			uh2=0ULL;
			for (p=1;p<N;p++)
			{
				uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
			
				fprintf(fid,"%hhu\n",uh+uh2+2*uinter);
			}
		}
	}
}
	

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_sq_c_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uh2,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
			utest = m ^ circ_bin_lshift(m,N,bin);
			uh2=0ULL;
			for (p=0;p<N;p++)
			{
				uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%hhu\n",uh+uh2+2*uinter);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_f_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uh2,uinter,xh,xh2;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
			utest = m ^ circ_bin_lshift(m,N,bin);
			uh2=0ULL;
			for (p=1;p<N;p++)
			{
				uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			xh2 = 0ULL;
			for (p=0;p<N;p++)
			{
				xh2 = xh2 + (((m & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%hhu %hhu\n",uh+uh2+2*uinter,xh+xh2);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_c_f_s_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
const unsigned char csize=CHAR_BIT;
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest;
unsigned char uh,uh2,uinter,xh,xh2;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_f_s_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
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
			utest = m ^ circ_bin_lshift(m,N,bin);
			uh2=0ULL;
			for (p=0;p<N;p++)
			{
				uh2 = uh2 + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			xh2 = 0ULL;
			for (p=0;p<N;p++)
			{
				xh2 = xh2 + (((m & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
			}
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				utest = n ^ m;
				uinter=0ULL;
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%hhu %hhu\n",uh+uh2+2*uinter,xh+xh2);
			}
		}
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}

