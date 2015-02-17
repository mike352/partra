#include "partra_genfuncs.h"
#include "partra_reductions.h"
#include "partra_ising.h"


/*****************************************************/
/**********Full Ising square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_f(const unsigned char N, char* dirname)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_sq_c_f(const unsigned char N, char* dirname)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_sq_f_f(const unsigned char N, char* dirname)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_sq_c_f(const unsigned char N, char* dirname)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*****************************************************/
/*******Full Ising triangular transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char i_tri_f_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_tri_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%d\n",(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_tri_c_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh;
char filename[256];

sprintf(filename,"%s/i_tri_c_%d.txt",dirname,N);
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
	    fprintf(fid,"%d\n",(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_f_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_tri_f_%d.txt",dirname,N);
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
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(m,N)))+uh),(xh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_c_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char uh,xh;
char filename[256];

sprintf(filename,"%s/if_tri_c_%d.txt",dirname,N);
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
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+bit_sum(n^circ_single_lshift(m,N))+uh),(xh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*****************************************************/
/********Reduced Ising square transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_r(const unsigned char N, char* dirname)
{
unsigned char umax=2*N-1;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_sq_f_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
				ctotal++;
				flip2 = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
				nn=n;
				for (p=0;p<flip+1;p++)
				{
					mm = m;
					for (q=0;q<flip2+1;q++)
					{
						melement[(bit_sum(nn^mm)+uh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu \n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char i_sq_c_r(const unsigned char N, char* dirname)
{
unsigned char umax=2*N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_sq_c_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum(n^circ_single_lshift(n,N));
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac3=lldiv(m,csize);
			if (((bitarray[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) //bitarray unique configuration
			{
				bitfrac4=lldiv(ctotal,csize);
				flip2 = (reflec[bitfrac4.quot]&(1<<bitfrac4.rem))>>bitfrac4.rem; //whether to reflect configuration
				nn=n;
				mm=m;
				for (p=0;p<flip+1;p++)
				{
					for (q=0;q<flip2+1;q++)
					{
						for (r=0;r<order[rtotal]+1;r++)
						{
							for (s=0;s<order[ctotal]+1;s++)
							{
								melement[(bit_sum(nn^mm)+uh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,melement[p]/((1ULL+flip2)*(order[ctotal]+1ULL)),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_sq_f_r(const unsigned char N, char* dirname)
{
unsigned char umax=2*N-1;
unsigned char xmax=N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((umax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (umax+1); n++)
{
	melement[n] = (unsigned char*) calloc(xmax+1, sizeof(unsigned char));
	if ((melement[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (m=0ULL;m<n;m++)
		{
			free((void*)melement[m]);
		}
		free((void*)melement);
		return 2;
	}
}


sprintf(filename,"%s/if_sq_f_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		xh = N-(bit_sum(n));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
				ctotal++;
				flip2 = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
				nn=n;
				for (p=0;p<flip+1;p++)
				{
					mm = m;
					for (q=0;q<flip2+1;q++)
					{
						melement[(bit_sum(nn^mm)+uh)][(xh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal,ctotal,melement[p][q]/(1ULL+flip2),p,q); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
for (n=0ULL;n<(umax+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_sq_c_r(const unsigned char N, char* dirname)
{
unsigned char umax=2*N;
unsigned char xmax=N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((umax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (umax+1); n++)
{
	melement[n] = (unsigned char*) calloc(xmax+1, sizeof(unsigned char));
	if ((melement[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (m=0ULL;m<n;m++)
		{
			free((void*)melement[m]);
		}
		free((void*)melement);
		return 2;
	}
}

sprintf(filename,"%s/if_sq_c_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum(n^circ_single_lshift(n,N));
		xh = N-(bit_sum(n));
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac3=lldiv(m,csize);
			if (((bitarray[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) //bitarray unique configuration
			{
				bitfrac4=lldiv(ctotal,csize);
				flip2 = (reflec[bitfrac4.quot]&(1<<bitfrac4.rem))>>bitfrac4.rem; //whether to reflect configuration
				nn=n;
				mm=m;
				for (p=0;p<flip+1;p++)
				{
					for (q=0;q<flip2+1;q++)
					{
						for (r=0;r<order[rtotal]+1;r++)
						{
							for (s=0;s<order[ctotal]+1;s++)
							{
								melement[((bit_sum(nn^mm))+uh)][(xh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,melement[p][q]/((1ULL+flip2)*(order[ctotal]+1ULL)),p,q); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
for (n=0ULL;n<(umax+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}



/*****************************************************/
/******Reduced Ising triangular transfer matrices*****/
/*****************************************************/

/*******************************/
unsigned char i_tri_f_r(const unsigned char N, char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
scanf("%s",option);
if (strcmp(option,"y")!=0)
{
	return 4;
}

unsigned char umax=3*N-2;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_tri_f_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
				ctotal++;
				flip2 = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
				nn=n;
				for (p=0;p<flip+1;p++)
				{
					mm = m;
					for (q=0;q<flip2+1;q++)
					{
						melement[(bit_sum(nn^mm)+bit_sum((~1ULL&nn)^(~1ULL&circ_single_lshift(mm,N)))+uh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu \n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char i_tri_c_r(const unsigned char N, char* dirname)
{
unsigned char umax=3*N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_tri_c_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum(n^circ_single_lshift(n,N));
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac3=lldiv(m,csize);
			if (((bitarray[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) //bitarray unique configuration
			{
				bitfrac4=lldiv(ctotal,csize);
				flip2 = (reflec[bitfrac4.quot]&(1<<bitfrac4.rem))>>bitfrac4.rem; //whether to reflect configuration
				nn=n;
				mm=m;
				for (p=0;p<flip+1;p++)
				{
					for (q=0;q<flip2+1;q++)
					{
						for (r=0;r<order[rtotal]+1;r++)
						{
							for (s=0;s<order[ctotal]+1;s++)
							{
								melement[(bit_sum(nn^mm)+bit_sum(nn^circ_single_lshift(mm,N))+uh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,melement[p]/((1ULL+flip2)*(order[ctotal]+1ULL)),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_tri_f_r(const unsigned char N, char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
scanf("%s",option);
if (strcmp(option,"y")!=0)
{
	return 4;
}

unsigned char umax=3*N-2;
unsigned char xmax=N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((umax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (umax+1); n++)
{
	melement[n] = (unsigned char*) calloc(xmax+1, sizeof(unsigned char));
	if ((melement[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (m=0ULL;m<n;m++)
		{
			free((void*)melement[m]);
		}
		free((void*)melement);
		return 2;
	}
}


sprintf(filename,"%s/if_tri_f_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		xh = N-(bit_sum(n));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
				ctotal++;
				flip2 = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
				nn=n;
				for (p=0;p<flip+1;p++)
				{
					mm = m;
					for (q=0;q<flip2+1;q++)
					{
						melement[(bit_sum(nn^mm)+bit_sum((~1ULL&nn)^(~1ULL&circ_single_lshift(mm,N)))+uh)][(xh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal,ctotal,melement[p][q]/(1ULL+flip2),p,q); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
for (n=0ULL;n<(umax+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_tri_c_r(const unsigned char N, char* dirname)
{
unsigned char umax=3*N;
unsigned char xmax=N;

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((umax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (umax+1); n++)
{
	melement[n] = (unsigned char*) calloc(xmax+1, sizeof(unsigned char));
	if ((melement[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (m=0ULL;m<n;m++)
		{
			free((void*)melement[m]);
		}
		free((void*)melement);
		return 2;
	}
}

sprintf(filename,"%s/if_tri_c_r_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum(n^circ_single_lshift(n,N));
		xh = N-(bit_sum(n));
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac3=lldiv(m,csize);
			if (((bitarray[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) //bitarray unique configuration
			{
				bitfrac4=lldiv(ctotal,csize);
				flip2 = (reflec[bitfrac4.quot]&(1<<bitfrac4.rem))>>bitfrac4.rem; //whether to reflect configuration
				nn=n;
				mm=m;
				for (p=0;p<flip+1;p++)
				{
					for (q=0;q<flip2+1;q++)
					{
						for (r=0;r<order[rtotal]+1;r++)
						{
							for (s=0;s<order[ctotal]+1;s++)
							{
								melement[(bit_sum(nn^mm)+bit_sum(nn^circ_single_lshift(mm,N))+uh)][(xh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,melement[p][q]/((1ULL+flip2)*(order[ctotal]+1ULL)),p,q); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
for (n=0ULL;n<(umax+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}