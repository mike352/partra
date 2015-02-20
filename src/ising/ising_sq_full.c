#include "partra_genfuncs.h"
#include "partra_ising.h"

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

/*****************************************************/
/**********Full Ising square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_sq_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%d\n",(bit_sum(n^m)+uh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_sq_c_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_sq_c_%d.txt",dirname,N);
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
unsigned char if_sq_f_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_sq_f_%d.txt",dirname,N);
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
unsigned char if_sq_c_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_sq_c_%d.txt",dirname,N);
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
