#include "partra_genfuncs.h"
#include "partra_reductions.h"
#include "partra_potts.h"


/*****************************************************/
/*****Reduced Potts triangular transfer matrices******/
/*****************************************************/

/*******************************/
unsigned char p2_tri_f_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
int scheck=scanf("%s",option);
if (scheck<=0)
{
	printf("ERROR: %s",strerror(errno));
	return 0;
}
if (strcmp(option,"y")!=0)
{
	return 4;
}
unsigned char umax = 3*N-2;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}


melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}


sprintf(filename,"%s/p_tri_f_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_bin_f(N,bin,&bitarray,&reflec,&total);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
						utest = nn ^ mm;
						uinter=0ULL;
						for (r=0;r<N;r++)
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						utest = nn ^ circ_bin_lshift(mm,N,bin);
						for (r=1;r<N;r++) //r starts at 1 for free b.c.
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[uh+uinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,p,melement[p]/(1ULL+flip2)); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}

/*******************************/
unsigned char p2_tri_c_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
unsigned char umax = 3*N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/p_tri_c_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_bin_c(N,bin,&bitarray,&reflec,&order,&total);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
								utest = nn ^ mm;
								uinter=0ULL;
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								utest = nn ^ circ_bin_lshift(mm,N,bin);
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[uh+uinter]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,p,melement[p]/((1ULL+flip2)*(order[ctotal]+1ULL))); //add normalization constant only to row vector
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

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;

}

/*******************************/
unsigned char pf2_tri_f_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
int scheck=scanf("%s",option);
if (scheck<=0)
{
	printf("ERROR: %s",strerror(errno));
	return 0;
}
if (strcmp(option,"y")!=0)
{
	return 4;
}
unsigned char umax = 3*N-2;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}


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


sprintf(filename,"%s/pf_tri_f_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_bin_f(N,bin,&bitarray,&reflec,&total);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh = 0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
						utest = nn ^ mm;
						uinter=0ULL;
						for (r=0;r<N;r++)
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						utest = nn ^ circ_bin_lshift(mm,N,bin);
						for (r=1;r<N;r++) //r starts at 1 for free b.c.
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[uh+uinter][xh]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<(xmax+1);q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal,ctotal,p,q,melement[p][q]/(1ULL+flip2)); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile %s created.",filename);
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
unsigned char pf2_tri_c_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
unsigned char umax = 3*N;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

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

sprintf(filename,"%s/pf_tri_c_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_simple_bin_c(N,bin,&bitarray,&reflec,&order,&total);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
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
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
								utest = nn ^ mm;
								uinter=0ULL;
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								utest = nn ^ circ_bin_lshift(mm,N,bin);
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[uh+uinter][xh]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<(xmax+1);q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,p,q,melement[p][q]/((1ULL+flip2)*(order[ctotal]+1ULL))); //add normalization constant only to row vector
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

printf("\nFile %s created.",filename);
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


/*******************************/
unsigned char p_tri_f_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
int scheck=scanf("%s",option);
if (scheck<=0)
{
	printf("ERROR: %s",strerror(errno));
	return 0;
}
if (strcmp(option,"y")!=0)
{
	return 4;
}
unsigned char umax = 3*N-2;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,nn,mm,sum;
unsigned char* pnums;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
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


sprintf(filename,"%s/p_tri_f_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_gen_bin_f(N,bin,&bitarray,&reflec,&total,pnums);
//printf("\nN=%llu: unique=%llu\n",N,total); 
free((void*)pnums);
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/


//Calculate transfer matrix
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
						utest = nn ^ mm;
						uinter=0ULL;
						for (r=0;r<N;r++)
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						utest = nn ^ circ_bin_lshift(mm,N,bin);
						for (r=1;r<N;r++) //r starts at 1 for free b.c.
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[uh+uinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,p,melement[p]/(1ULL+flip2)); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char p_tri_c_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
unsigned char umax = 3*N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
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

sprintf(filename,"%s/p_tri_c_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_gen_bin_c(N,bin,&bitarray,&reflec,&order,&total,pnums);
//printf("\nN=%llu: unique=%llu\n",N,total); 
free((void*)pnums);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
								utest = nn ^ mm;
								uinter=0ULL;
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								utest = nn ^ circ_bin_lshift(mm,N,bin);
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[uh+uinter]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,p,melement[p]/((1ULL+flip2)*(order[ctotal]+1ULL))); //add normalization constant only to row vector
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

printf("\nFile %s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char pf_tri_f_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
char option[256];
printf("\nThe full transfer matrix for free row boundary conditions does not have a direct sum in terms of parity sectors. Therefore the reduced transfer matrix is not a valid reflection symmetric sector. Continue anyway? (y,n): ");
int scheck=scanf("%s",option);
if (scheck<=0)
{
	printf("ERROR: %s",strerror(errno));
	return 0;
}
if (strcmp(option,"y")!=0)
{
	return 4;
}

unsigned char umax = 3*N-2;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,nn,mm,sum;
unsigned char* pnums;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

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


sprintf(filename,"%s/pf_tri_f_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_gen_bin_f(N,bin,&bitarray,&reflec,&total,pnums);
//printf("\nN=%llu: unique=%llu\n",N,total); 
free((void*)pnums);
if (flag!=0)
{
	fclose(fid);
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/


//Calculate transfer matrix
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		rtotal++;
		utest = n ^ circ_bin_lshift(n,N,bin);
		uh=0ULL;
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uh = uh + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		xh = 0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
						utest = nn ^ mm;
						uinter=0ULL;
						for (r=0;r<N;r++)
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						utest = nn ^ circ_bin_lshift(mm,N,bin);
						for (r=1;r<N;r++) //r starts at 1 for free b.c.
						{
							uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[uh+uinter][xh]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<(xmax+1);q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal,ctotal,p,q,melement[p][q]/(1ULL+flip2)); //add normalization constant only to row vector
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile %s created.",filename);
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
unsigned char pf_tri_c_r_file(const unsigned char N, const unsigned long long Q, const char* dirname)
{
unsigned char umax = 3*N;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long utest,rtotal=0ULL,ctotal=0ULL;
unsigned char uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

pnums = (unsigned char*) calloc((1ULL<<bin*N)/csize+1ULL,sizeof(unsigned char));
if ((pnums==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

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

sprintf(filename,"%s/pf_tri_c_r_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = red_gen_bin_c(N,bin,&bitarray,&reflec,&order,&total,pnums);
//printf("\nN=%llu: unique=%llu\n",N,total); 
free((void*)pnums);
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
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
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
		bitfrac2=lldiv(rtotal,csize);
		flip = (reflec[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N*bin);m++)
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
								utest = nn ^ mm;
								uinter=0ULL;
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								utest = nn ^ circ_bin_lshift(mm,N,bin);
								for (t=0;t<N;t++)
								{
									uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[uh+uinter][xh]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<(xmax+1);q++)
					{
						if (melement[p][q]!=0)
						{
							fprintf(fid,"%llu %llu %llu %llu %llu\n",rtotal+1ULL,ctotal+1ULL,p,q,melement[p][q]/((1ULL+flip2)*(order[ctotal]+1ULL))); //add normalization constant only to row vector
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

printf("\nFile %s created.",filename);
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