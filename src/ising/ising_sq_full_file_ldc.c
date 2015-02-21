#include "partra_genfuncs.h"
#include "partra_ising_ldc.h"

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

/*****************************************************/
/**********Full Ising square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_f_file_ldc(const unsigned char N, const char* dirname,long double complex uc)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_sq_f_f_%d_%#f_I%#f.txt",dirname,N,creal(uc),cimag(uc));
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%#.32f %#.32f\n",creal(cpowl(uc,bit_sum(n^m)+uh)),cimag(cpowl(uc,bit_sum(n^m)+uh)));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_sq_c_f_file_ldc(const unsigned char N, const char* dirname,long double complex uc)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_sq_c_f_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum(n^circ_single_lshift(n,N));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d\n",(bit_sum(n^m)+uh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_sq_f_f_file_ldc(const unsigned char N, const char* dirname, const char* which, ...)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

unsigned char u=0,x=0;
long double complex uc=0,xc=0;
va_list vl;
va_start(vl,which);
if (strcmp(which,"u")==0)
{
	u=1;
	uc=va_arg(vl,long double complex);
}
else if (strcmp(which,"x")==0)
{
	x=1;
	xc=va_arg(vl,long double complex);
}
else if ((strcmp(which,"ux")==0))
{
	u=1;x=1;
	uc=va_arg(vl,long double complex);
	xc=va_arg(vl,long double complex);
}
else if ((strcmp(which,"xu")==0))
{
	u=1;x=1;
	xc=va_arg(vl,long double complex);
	uc=va_arg(vl,long double complex);
}
else
{
	printf("Wrong numerical substitution selection. Use only \"u\", \"x\", or \"ux\" (or \"xu\").\n");
	return 3;
}
va_end(vl);

sprintf(filename,"%s/if_sq_f_f_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	xh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+uh),(xh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_sq_c_f_file_ldc(const unsigned char N, const char* dirname, const char* which, ...)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

unsigned char u=0,x=0;
long double complex uc=0,xc=0;
va_list vl;
va_start(vl,which);
if (strcmp(which,"u")==0)
{
	u=1;
	uc=va_arg(vl,long double complex);
}
else if (strcmp(which,"x")==0)
{
	x=1;
	xc=va_arg(vl,long double complex);
}
else if ((strcmp(which,"ux")==0))
{
	u=1;x=1;
	uc=va_arg(vl,long double complex);
	xc=va_arg(vl,long double complex);
}
else if ((strcmp(which,"xu")==0))
{
	u=1;x=1;
	xc=va_arg(vl,long double complex);
	uc=va_arg(vl,long double complex);
}
else
{
	printf("Wrong numerical substitution selection. Use only \"u\", \"x\", or \"ux\" (or \"xu\").\n");
	return 3;
}
va_end(vl);

sprintf(filename,"%s/if_sq_c_f_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	uh = bit_sum(n^circ_single_lshift(n,N));
	xh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+uh),(xh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}
