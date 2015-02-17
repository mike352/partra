#include "partra_genfuncs.h"
#include "partra_reductions.h"

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

/*******************************/
unsigned char red_simple_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
{
unsigned long long n,m;
const unsigned char csize=CHAR_BIT;
lldiv_t bitfrac, bitfrac2;	

*bitarray = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}

//Find unique numbers, and whether reflection symmetric
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		m = bit_reflection(n,N);
		if (m!=n)
		{
			bitfrac2=lldiv(m,csize);
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
			(*reflec)[bitfrac.quot]=(*reflec)[bitfrac.quot] + (1<<bitfrac.rem);
		}
		(*total)++;
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec's bits which are 1 represent the non-reflection symmetric unique configurations

return 0;
}


/*******************************/
unsigned char red_simple_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=CHAR_BIT;

*bitarray = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

//Using A000029 sequence to minimize size of array
if (N<=10)
{
	*reflec = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
}
else if (N<18)
{
	*reflec = (unsigned char*) calloc((1ULL<<(N-4))/csize+1ULL,sizeof(unsigned char));
}
else if (N<33)
{
	*reflec = (unsigned char*) calloc((1ULL<<(N-5))/csize+1ULL,sizeof(unsigned char));
}
else
{
	*reflec = (unsigned char*) calloc((1ULL<<(N-6))/csize+1ULL,sizeof(unsigned char));
}
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}
	
//Using A000029 sequence to minimize size of array
if (N<=10)
{	
	*order = (unsigned char*) calloc(1ULL<<N,sizeof(unsigned char));
}
else if (N<18)
{
	*order = (unsigned char*) calloc(1ULL<<(N-4),sizeof(unsigned char));
}
else if (N<33)
{
	*order = (unsigned char*) calloc(1ULL<<(N-5),sizeof(unsigned char));
}
else
{
	*order = (unsigned char*) calloc(1ULL<<(N-6),sizeof(unsigned char));
}
if ((*order==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}


//Find unique numbers, their order, and whether reflection symmetric
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		m = circ_single_lshift(n,N);
		while (m!=n)
		{
			(*order)[*total]++; //order = 0 if unique. order=1 if two exist, etc.
			bitfrac2=lldiv(m,csize); 
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem); 
			m = circ_single_lshift(m,N); 
		}
		(*total)++; 
	}
	m=bit_reflection(n,N); 
	bitfrac3=lldiv(m,csize);
	if (((((*bitarray)[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) & (m>n))
	{
		(*bitarray)[bitfrac3.quot]=(*bitarray)[bitfrac3.quot] + (1<<bitfrac3.rem);
		bitfrac4=lldiv(*total-1ULL,csize);
		(*reflec)[bitfrac4.quot]=(*reflec)[bitfrac4.quot] + (1<<bitfrac4.rem);
		p = circ_single_lshift(m,N);
		while (p!=m)
		{
			bitfrac2=lldiv(p,csize);
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
			p = circ_single_lshift(p,N);
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec is a bitfield but of size total, indexed by unique configurations, which are 1 represent the on-reflection symmetric configurations
//unique.order is an array of size total and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.


//Save memory by reducing the size of the reflec array
*reflec = (unsigned char*) realloc(*reflec,((*total)/csize+1ULL)*sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*reflec);
	free((void*)*bitarray);
	free((void*)*order);
	return 2;
}

//Save memory by reducing the size of the order array
*order = (unsigned char*) realloc(*order,(*total)*sizeof(unsigned char));
if ((order==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*order);
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}

return 0;
}


/*******************************/
unsigned char red_simple_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
{
unsigned long long n,m;
const unsigned char csize=CHAR_BIT;
lldiv_t bitfrac, bitfrac2;	

*bitarray = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}

//Find unique numbers, and whether reflection symmetric
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		m = bit_reflection_bin(n,N,bin);
		if (m!=n)
		{
			bitfrac2=lldiv(m,csize);
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
			(*reflec)[bitfrac.quot]=(*reflec)[bitfrac.quot] + (1<<bitfrac.rem);
		}
		(*total)++;
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec's bits which are 1 represent the non-reflection symmetric unique configurations

return 0;
}


/*******************************/
unsigned char red_simple_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=CHAR_BIT;

*bitarray = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}
	
*order = (unsigned char*) calloc(1ULL<<N*bin,sizeof(unsigned char));
if ((*order==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}


//Find unique numbers, their order, and whether reflection symmetric
for (n=0;n<(1ULL<<(N*bin));n++)
{
	bitfrac=lldiv(n,csize);
	if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		m = circ_bin_lshift(n,N,bin);
		while (m!=n)
		{//printf("By\n");
			(*order)[*total]++; //order = 0 if unique. order=1 if two exist, etc.
			bitfrac2=lldiv(m,csize); 
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
			m = circ_bin_lshift(m,N,bin);
		}
		(*total)++; 
	}
	m=bit_reflection_bin(n,N,bin); 
	bitfrac3=lldiv(m,csize);
	if (((((*bitarray)[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) & (m>n))
	{
		(*bitarray)[bitfrac3.quot]=(*bitarray)[bitfrac3.quot] + (1<<bitfrac3.rem);
		bitfrac4=lldiv(*total-1ULL,csize);
		(*reflec)[bitfrac4.quot]=(*reflec)[bitfrac4.quot] + (1<<bitfrac4.rem);
		p = circ_bin_lshift(m,N,bin);
		while (p!=m)
		{
			bitfrac2=lldiv(p,csize);
			(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
			p = circ_bin_lshift(p,N,bin);
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec is a bitfield but of size total, indexed by unique configurations, which are 1 represent the on-reflection symmetric configurations
//unique.order is an array of size total and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.


//Save memory by reducing the size of the reflec array
*reflec = (unsigned char*) realloc(*reflec,((*total)/csize+1ULL)*sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*reflec);
	free((void*)*bitarray);
	free((void*)*order);
	return 2;
}

//Save memory by reducing the size of the order array
*order = (unsigned char*) realloc(*order,(*total)*sizeof(unsigned char));
if ((order==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*order);
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}

return 0;
}


/*******************************/
unsigned char red_gen_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m;
const unsigned char csize=CHAR_BIT;
lldiv_t bitfrac, bitfrac2;	

*bitarray = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}

//Find unique numbers, and whether reflection symmetric
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
		{
			m = bit_reflection(n,N);
			if (m!=n)
			{
				bitfrac2=lldiv(m,csize);
				if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
				{
					(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
					(*reflec)[bitfrac.quot]=(*reflec)[bitfrac.quot] + (1<<bitfrac.rem);
				}
				else
				{
					printf("\nERROR: Cannot reduce transfer matrix because reflected configurations are not valid configurations.");
					free((void*)*bitarray);
					free((void*)*reflec);
					return 4;
				}
			}
			(*total)++;
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec's bits which are 1 represent the non-reflection symmetric unique configurations

//Remove invalid configurations from bitarray
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==1)
	{
		(*bitarray)[bitfrac.quot]=(*bitarray)[bitfrac.quot] + (1<<bitfrac.rem);
	}
}

return 0;
}


/*******************************/
unsigned char red_gen_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=CHAR_BIT;

*bitarray = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}
	
*order = (unsigned char*) calloc(1ULL<<N,sizeof(unsigned char));
if ((*order==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}


//Find unique numbers, their order, and whether reflection symmetric
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
		{
			m = circ_single_lshift(n,N);
			while (m!=n)
			{
				bitfrac2=lldiv(m,csize); 
				if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
				{
					(*order)[*total]++; //order = 0 if unique. order=1 if two exist, etc.
					(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem); 
					m = circ_single_lshift(m,N); 
				}
				else
				{
					printf("\nERROR: Cannot reduce transfer matrix because rotated configurations are not valid configurations.");
					free((void*)*bitarray);
					free((void*)*reflec);
					free((void*)*order);
					return 4;
				}
			}
			(*total)++; 
		}
		m=bit_reflection(n,N); 
		bitfrac3=lldiv(m,csize);
		if (((numbers[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0)
		{
			if (((((*bitarray)[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) & (m>n))
			{
				(*bitarray)[bitfrac3.quot]=(*bitarray)[bitfrac3.quot] + (1<<bitfrac3.rem);
				bitfrac4=lldiv(*total-1ULL,csize);
				(*reflec)[bitfrac4.quot]=(*reflec)[bitfrac4.quot] + (1<<bitfrac4.rem);
				p = circ_single_lshift(m,N);
				while (p!=m)
				{
					bitfrac2=lldiv(p,csize);
					if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
					{
						(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
						p = circ_single_lshift(p,N);
					}
					else
					{
						printf("\nERROR: Cannot reduce transfer matrix because reflected and rotated configurations are not valid configurations.");
						free((void*)*bitarray);
						free((void*)*reflec);
						free((void*)*order);
						return 4;
					}
				}
			}
		}
		else
		{
			printf("\nERROR: Cannot reduce transfer matrix because reflected configurations are not valid configurations.");
			free((void*)*bitarray);
			free((void*)*reflec);
			free((void*)*order);
			return 4;
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec is a bitfield but of size total, indexed by unique configurations, which are 1 represent the on-reflection symmetric configurations
//unique.order is an array of size total and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.

//Remove invalid configurations from bitarray
for (n=0;n<(1ULL<<N);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==1)
	{
		(*bitarray)[bitfrac.quot]=(*bitarray)[bitfrac.quot] + (1<<bitfrac.rem);
	}
}

//Save memory by reducing the size of the reflec array
*reflec = (unsigned char*) realloc(*reflec,((*total)/csize+1ULL)*sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*reflec);
	free((void*)*bitarray);
	free((void*)*order);
	return 2;
}

//Save memory by reducing the size of the order array
*order = (unsigned char*) realloc(*order,(*total)*sizeof(unsigned char));
if ((order==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*order);
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}

return 0;
}


/*******************************/
unsigned char red_gen_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m;
const unsigned char csize=CHAR_BIT;
lldiv_t bitfrac, bitfrac2;	

*bitarray = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}

//Find unique numbers, and whether reflection symmetric
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
		{
			m = bit_reflection_bin(n,N,bin);
			if (m!=n)
			{
				bitfrac2=lldiv(m,csize);
				if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
				{
					(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
					(*reflec)[bitfrac.quot]=(*reflec)[bitfrac.quot] + (1<<bitfrac.rem);
				}
				else
				{
					printf("\nERROR: Cannot reduce transfer matrix because reflected configurations are not valid configurations.");
					free((void*)*bitarray);
					free((void*)*reflec);
					return 4;
				}
			}
			(*total)++;
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec's bits which are 1 represent the non-reflection symmetric unique configurations


//Remove invalid configurations from bitarray
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==1)
	{
		(*bitarray)[bitfrac.quot]=(*bitarray)[bitfrac.quot] + (1<<bitfrac.rem);
	}
}

return 0;
}


/*******************************/
unsigned char red_gen_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=CHAR_BIT;

*bitarray = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*bitarray==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

*reflec = (unsigned char*) calloc((1ULL<<N*bin)/csize+1ULL,sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	return 2;
}
	
*order = (unsigned char*) calloc(1ULL<<N*bin,sizeof(unsigned char));
if ((*order==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}


//Find unique numbers, their order, and whether reflection symmetric
for (n=0;n<(1ULL<<(N*bin));n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
	{
		if ((((*bitarray)[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==0)
		{
			m = circ_bin_lshift(n,N,bin);
			while (m!=n)
			{
				bitfrac2=lldiv(m,csize); 
				if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
				{
					(*order)[*total]++; //order = 0 if unique. order=1 if two exist, etc.
					(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
					m = circ_bin_lshift(m,N,bin);
				}
				else
				{
					printf("\nERROR: Cannot reduce transfer matrix because rotated configurations are not valid configurations.");
					free((void*)*bitarray);
					free((void*)*reflec);
					free((void*)*order);
					return 4;
				}
			}
			(*total)++; 
		}
		m=bit_reflection_bin(n,N,bin); 
		bitfrac3=lldiv(m,csize);
		if (((numbers[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0)
		{
			if (((((*bitarray)[bitfrac3.quot]&(1<<bitfrac3.rem))>>bitfrac3.rem)==0) & (m>n))
			{
				(*bitarray)[bitfrac3.quot]=(*bitarray)[bitfrac3.quot] + (1<<bitfrac3.rem);
				bitfrac4=lldiv(*total-1ULL,csize);
				(*reflec)[bitfrac4.quot]=(*reflec)[bitfrac4.quot] + (1<<bitfrac4.rem);
				p = circ_bin_lshift(m,N,bin);
				while (p!=m)
				{
					bitfrac2=lldiv(p,csize);
					if (((numbers[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
					{
						(*bitarray)[bitfrac2.quot]=(*bitarray)[bitfrac2.quot] + (1<<bitfrac2.rem);
						p = circ_bin_lshift(p,N,bin);
					}
					else
					{
						printf("\nERROR: Cannot reduce transfer matrix because reflected and rotated configurations are not valid configurations.");
						free((void*)*bitarray);
						free((void*)*reflec);
						free((void*)*order);
						return 4;
					}
				}
			}
		}
		else
		{
			printf("\nERROR: Cannot reduce transfer matrix because reflected configurations are not valid configurations.");
			free((void*)*bitarray);
			free((void*)*reflec);
			free((void*)*order);
			return 4;
		}
	}
}
//unique.bitarray's zero bits represent the unique configurations
//unique.reflec is a bitfield but of size total, indexed by unique configurations, which are 1 represent the on-reflection symmetric configurations
//unique.order is an array of size total and gives order of rotations. order=0 for no extra rotations, order=1 for one rotation, etc.


//Remove invalid configurations from bitarray
for (n=0;n<(1ULL<<N*bin);n++)
{
	bitfrac=lldiv(n,csize);
	if (((numbers[bitfrac.quot]&(1<<bitfrac.rem))>>bitfrac.rem)==1)
	{
		(*bitarray)[bitfrac.quot]=(*bitarray)[bitfrac.quot] + (1<<bitfrac.rem);
	}
}

//Save memory by reducing the size of the reflec array
*reflec = (unsigned char*) realloc(*reflec,((*total)/csize+1ULL)*sizeof(unsigned char));
if ((*reflec==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*reflec);
	free((void*)*bitarray);
	free((void*)*order);
	return 2;
}

//Save memory by reducing the size of the order array
*order = (unsigned char*) realloc(*order,(*total)*sizeof(unsigned char));
if ((order==NULL))
{
	printf("\nERROR: Could not re-allocate memory.");
	free((void*)*order);
	free((void*)*bitarray);
	free((void*)*reflec);
	return 2;
}

return 0;
}