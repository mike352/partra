#include "partra_genfuncs.h"
#include "partra_reductions.h"
#include "partra_ising.h"

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

/*****************************************************/
/********Reduced Ising square transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char i_sq_f_r(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
unsigned char umax=2*N-1;
sprintf(filename,"i_sq_f_r_%d.txt",N);

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;
unsigned char count;

melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

msize[0] = total;
flag = matrix_alloc(matrix,msize,N);
if (flag!=0)
{
	free((void*)bitarray);
	free((void*)reflec);
	free((void*)melement);
	return flag;
}

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
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
				count=0ULL;
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						//add normalization constant only to row vector
						(*matrix)[rtotal][ctotal][1][msize[1]*count]=p;
						(*matrix)[rtotal][ctotal][1][msize[1]*count+1]=melement[p]/(1ULL+flip2);
						count++;
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
				(*matrix)[rtotal][ctotal][0][0]=count;
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char i_sq_c_r(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=2ULL;
unsigned char umax=2*N;
sprintf(filename,"i_sq_c_r_%d.txt",N);

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,rtotal=0ULL,ctotal=0ULL;
unsigned char count;

melement = (unsigned char*) calloc((umax+1),sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

msize[0] = total;
flag = matrix_alloc(matrix,msize,N);
if (flag!=0)
{
	free((void*)bitarray);
	free((void*)reflec);
	free((void*)melement);
	return flag;
}

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
				count=0ULL;
				for (p=0;p<umax+1;p++)
				{
					if (melement[p]!=0)
					{
						//add normalization constant only to row vector
						(*matrix)[rtotal][ctotal][1][msize[1]*count]=p;
						(*matrix)[rtotal][ctotal][1][msize[1]*count+1]=melement[p]/((1ULL+flip2)*(order[ctotal]+1ULL));
						count++;
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
				(*matrix)[rtotal][ctotal][0][0]=count;
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_sq_f_r(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
unsigned char umax=2*N-1;
unsigned char xmax=N;
sprintf(filename,"if_sq_f_r_%d.txt",N);

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;
unsigned char count;

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


//Compute unique configurations
flag = red_simple_f(N,&bitarray,&reflec,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	return flag;
}
/*Conventions: 
bitarray's zero bits are unique configurations
reflec's 1 bits are not reflection symmetric*/

msize[0] = total;
flag = matrix_alloc(matrix,msize,N);
if (flag!=0)
{
	free((void*)bitarray);
	free((void*)reflec);
	free((void*)melement);
	return flag;
}

//Calculate transfer matrix
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((bitarray[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0) //bitarray unique configuration
	{
		uh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		xh = N-(bit_sum(n));
		flip = (reflec[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem; //whether to reflect configuration
		for (m=0;m<(1ULL<<N);m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((bitarray[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0) //bitarray unique configuration
			{
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
				count=0ULL;
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							//add normalization constant only to row vector
							(*matrix)[rtotal][ctotal][1][msize[1]*count]=p;
							(*matrix)[rtotal][ctotal][1][msize[1]*count+1]=q;
							(*matrix)[rtotal][ctotal][1][msize[1]*count+2]=melement[p][q]/(1ULL+flip2); 
							count++;
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
						}
					}
				}
				(*matrix)[rtotal][ctotal][0][0]=count;
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

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
unsigned char if_sq_c_r(unsigned char***** matrix, unsigned long long* msize, char* filename, const unsigned char N)
{
msize[1]=3ULL;
unsigned char umax=2*N;
unsigned char xmax=N;
sprintf(filename,"if_sq_c_r_%d.txt",N);

unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
const unsigned char csize=CHAR_BIT;
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,xh,rtotal=0ULL,ctotal=0ULL;
unsigned char count;

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

//Compute unique configurations
flag = red_simple_c(N,&bitarray,&reflec,&order,&total);
//printf("\nN=%llu: unique=%llu\n",N,total); 
if (flag!=0)
{
	return flag;
}
/*Conventions: 
bitarray is a bitfield of size 2^N, zero bits are unique configurations
reflec is a bitfield of size total, indexed by unique configurations, 1 bits are not reflection symmetric
order is an array of size total indexed by unique configuration and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.*/

msize[0] = total;
flag = matrix_alloc(matrix,msize,N);
if (flag!=0)
{
	free((void*)bitarray);
	free((void*)reflec);
	free((void*)melement);
	return flag;
}

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
				count=0ULL;
				for (p=0;p<umax+1;p++)
				{
					for (q=0;q<xmax+1;q++)
					{
						if (melement[p][q]!=0)
						{
							//add normalization constant only to row vector
							(*matrix)[rtotal][ctotal][1][msize[1]*count]=p;
							(*matrix)[rtotal][ctotal][1][msize[1]*count+1]=q;
							(*matrix)[rtotal][ctotal][1][msize[1]*count+2]=melement[p][q]/((1ULL+flip2)*(order[ctotal]+1ULL)); 
							melement[p][q]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
							count++;
						}
					}
				}
				(*matrix)[rtotal][ctotal][0][0]=count;
				ctotal++;
			}
		}
		ctotal=0ULL;
		rtotal++;
	}
}

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
