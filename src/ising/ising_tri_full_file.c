#include "partra_genfuncs.h"
#include "partra_ising.h"


/*****************************************************/
/*******Full Ising triangular transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char i_tri_f_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_tri_f_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%hhu\n",(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_tri_c_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_tri_c_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%hhu\n",(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh));
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_f_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_tri_f_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%hhu %hhu\n",(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh),xh);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_c_f_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_tri_c_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%hhu %hhu\n",(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh),xh);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}



/*****************************************************/
/******************SYMMETRIC MATRICES*****************/
/*****************************************************/

/*******************************/
unsigned char i_tri_f_f_s_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,uh2;
char filename[256];

sprintf(filename,"%s/i_tri_f_f_s_%d.txt",dirname,N);
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
	    uh2 = bit_sum((~1ULL&m)^(~1ULL&circ_single_lshift(m,N)));
	    fprintf(fid,"%hhu\n",2*(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N))))+uh+uh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_tri_c_f_s_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,uh2;
char filename[256];

sprintf(filename,"%s/i_tri_c_f_s_%d.txt",dirname,N);
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
	    uh2 = bit_sum(m^circ_single_lshift(m,N));
	    fprintf(fid,"%hhu\n",2*(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N)))+uh+uh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_f_f_s_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,uh2,xh,xh2;
char filename[256];

sprintf(filename,"%s/if_tri_f_f_s_%d.txt",dirname,N);
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
	    uh2 = bit_sum((~1ULL&m)^(~1ULL&circ_single_lshift(m,N)));
	    xh2 = N-(bit_sum(m));
	    fprintf(fid,"%hhu %hhu\n",2*(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N))))+uh+uh2,xh+xh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_c_f_s_file(const unsigned char N, const char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,uh2,xh,xh2;
char filename[256];

sprintf(filename,"%s/if_tri_c_f_s_%d.txt",dirname,N);
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
	    uh2 = bit_sum(m^circ_single_lshift(m,N));
	    xh2 = N-(bit_sum(m));
	    fprintf(fid,"%hhu %hhu\n",2*(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N)))+uh+uh2,xh+xh2);
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
return 0;
}