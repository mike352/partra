/* Options:
1. Ising or Ising in a field or Potts or Potts in a field
2. row boundary conditions 
3. full or reduced transfer matrix
4. row size

Output file:
All data is output to a file in a (constructed) directory "data". Filename is displayed on screen after successful completion of calculations and writing of data to file. 

Output file format legend: 
r = row
c = column
n: x^n = exp(-2E/kT)^n, energy exponent
m: u^m = exp(-2H/kT)^m, field exponent. In q-Potts, H-field only acts on q=1
A: A x^n u^m, constant 

Output file format:
1. Full transfer matrix:
n m
(r and c are assumed)
When H=0, m is not output

2. Reduced transfer matrix
r c A n m
When H=0, m is not output

To Add:
1. Change bounds in for loops for reduced tm functions
2. Change output filenames
3. Change function names to specify lattice type
4. Add triangular lattice versions
5. Add api versions
6. Write a portable bit architecture tester 
7. Write a test suite


*/

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
//#include <direct.h> //uncomment for Windows
#include <sys/stat.h> //uncomment for Linux



/*This function performs a circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_single_lshift(unsigned long long number, const unsigned char width)
{
	//(shift left of number without Nth digit) + (Nth digit shifted to the front)
	number = ((~(~0ULL<<(width-1))& number)<<1) + (((1ULL<<(width-1))& number)>>(width-1));
	return number;
}

/*This function performs a binned circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_bin_lshift(unsigned long long number, const unsigned char N, const unsigned char bin)
{
	unsigned char n;
	for (n=0;n<bin;n++)
	{
		//(shift left of number without width digit) + (width digit shifted to the front)
		number = ((~(~0ULL<<(N*bin-1))& number)<<1) + (((1ULL<<(N*bin-1))& number)>>(N*bin-1));
	}
	return number;
}

/*This function reflects all the bits for a number of a particular width*/
unsigned long long bit_reflection(unsigned long long number, const unsigned char N)
{
	unsigned long long reflection=number;
	unsigned char s = N-1;
	//Bit twiddling from https://graphics.stanford.edu/~seander/bithacks.html
	for (number >>= 1; number; number >>= 1)
	{   
		reflection <<= 1;
		reflection |= number & 1ULL;
		s--;
	}
	reflection <<= s; 
	return (reflection & ~((~0ULL<<(N-1))<<1)); //added to algorithm to clear any set bits beyond N, done in two steps to handle case when N equals the size of long long
}

/*This function reflects all the bits for a number of a particular width and bin size bin*/
unsigned long long bit_reflection_bin(unsigned long long number, const unsigned char N, const unsigned char bin)
{
	unsigned char ii;
	unsigned long long reflection=0ULL;
	for (ii=0;ii<N;ii++)
	{
		//reflection=reflection + (select out ii-th bin and shift to starting position, shift appropriately)
		reflection = reflection + (((number & (((1ULL<<bin)-1ULL)<<bin*ii))>>bin*ii) << (bin*(N-1-ii)));
	}
	return reflection;
}

/*This function adds up all of the set bits in a number x*/
unsigned char bit_sum(unsigned long long x)
{
	unsigned char sum;
	//Standard bit twiddling
	for (sum=0;x;sum++) x &= x-1ULL;
	return sum;
}


//Main function declarations
unsigned char i_f(const unsigned char, char*); //Ising full transfer matrix, free row b.c.
unsigned char i_c(const unsigned char, char*); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_f(const unsigned char, char*); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_c(const unsigned char, char*); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_f_p(const unsigned char, char*); //Ising reduced transfer matrix, free row b.c.
unsigned char i_c_0p(const unsigned char, char*); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_f_p(const unsigned char, char*); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_c_0p(const unsigned char, char*); //Ising in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c.
unsigned char p_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_f_p(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_c_0p(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_f_p(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_c_0p(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p2_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_f_p(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_c_0p(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_f_p(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_c_0p(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

unsigned char simple_red_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N
unsigned char simple_red_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N
unsigned char simple_red_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N*bin
unsigned char simple_red_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N*bin

unsigned char arb_red_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char arb_red_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char arb_red_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array
unsigned char arb_red_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

int main(void)
{
unsigned char row_max_size=0;
unsigned long long test;
unsigned char N;
char Qs[256];
unsigned long long Q=0ULL;
unsigned char bin=0;
char dirname[256];
int dcheck; //Needs to be of type int for use with mkdir
time_t tic;
time_t toc;
double totaltime;
char option1[3];
char option2[2];
char option3[2];
unsigned char flag;


//Discover max row size
test = 1ULL;
while (test!=0ULL)
{
	row_max_size++;
	test = test << 1;
}

//Create directory. I haven't tested this in Windows.
sprintf(dirname,"data");
dcheck = mkdir(dirname,S_IRWXU | S_IRWXG | S_IRWXO);
if ((dcheck==-1)&(errno!=EEXIST))
{
	printf("\nERROR: Could not create output directory. %s\n",strerror(errno));
	return 0;
}

//Ising or Ising in a field
printf("\nIsing or Ising in a field or Potts or Potts in a field  (i,if,p,pf): ");
scanf("%s",option1);
if ((strcmp(option1,"i")!=0)&(strcmp(option1,"if")!=0)&(strcmp(option1,"p")!=0)&(strcmp(option1,"pf")!=0))
{
	printf("\nERROR: Wrong input.");
	return 0;
}

//Input q value, checking carefully, since very large numbers are allowed. 
if ((strcmp(option1,"p")==0)|(strcmp(option1,"pf")==0))
{
	printf("Potts q value  (3 to 2^%d-1): ",row_max_size);
	scanf("%255s",Qs);
	if (strspn(Qs,"1234567890")<strlen(Qs))
	{
		printf("\nERROR: Invalid q value.\n"); //prevents negative values which would get reinterpreted, floats, and others
		return 0;
	}
	Q=strtoll(Qs,NULL,10);
	if ((Q==LONG_MAX)|(Q==LONG_MIN))
	{
		printf("\nERROR: Given q value is too large. %s.\n",strerror(errno));
		return 0;
	}
	else if (Q==0)
	{
		printf("\nERROR: Invalid q value.\n"); //problem converting string Qs to unsigned long long Q
		return 0;
	}
	else if (Q<3ULL)
	{
		printf("\nERROR: q should be greater than 2.\n");
		return 0;
	}
	else if (Q>=3ULL)
	{
		while((1ULL<<bin)<Q)
		{
			bin++; //calculate bit array bin size
		}
	}
}

//Ask for row boundary condition
printf("Row boundary condition (f,c): ");
scanf("%s",option2);
if ((strcmp(option2,"f")!=0)&(strcmp(option2,"c")!=0))
{
	printf("\nERROR: Wrong input.");
	return 0;
}

//Ask for transfer matrix reduction
printf("Full or reduced transfer matrix (f,r): ");
scanf("%s",option3);
if ((strcmp(option3,"f")!=0)&(strcmp(option3,"r")!=0))
{
	printf("\nERROR: Wrong input.");
	return 0;
}


//Ask for row size
if ((strcmp(option1,"i")==0)|(strcmp(option1,"if")==0))
{
	printf("Row size (1 to %hhu): ",row_max_size); 
	scanf("%hhu",&N);
	if (N<1)
	{
		printf("\nERROR: Row size should be greater than 0.\n");
		return 0;
	}
	else if (N>row_max_size)
	{
		printf("\nERROR: Your machine can only do %d-bit calcuations.\n       Limit row size to %d.\n",row_max_size,row_max_size);
		return 0;
	}
}
else if ((strcmp(option1,"p")==0)|(strcmp(option1,"pf")==0))
{
	printf("Row size (1 to %hhu): ",row_max_size/bin); 
	scanf("%hhu",&N);
	if (N<1)
	{
		printf("\nERROR: Row size should be greater than 0.\n");
		return 0;
	}
	else if (N>row_max_size/bin)
	{
		printf("\nERROR: Your machine can only do %d-bit calcuations.\n       Limit row size to %d.\n",row_max_size,row_max_size/(bin*N));
		return 0;
	}
}

//Choose a function
time(&tic);
if (strcmp(option1,"i")==0)
{
	if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
	{
		flag = i_f(N,dirname);
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
	{
		flag = i_c(N,dirname);
	}
	else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
	{
		flag = i_f_p(N,dirname);
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
	{
		flag = i_c_0p(N,dirname);
	}
}
else if (strcmp(option1,"if")==0)
{
	if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
	{
		flag = if_f(N,dirname);
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
	{
		flag = if_c(N,dirname);
	}
	else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
	{
		flag = if_f_p(N,dirname);
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
	{
		flag = if_c_0p(N,dirname);
	}
}
else if (strcmp(option1,"p")==0)
{
	if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
	{
		if (Q==(1<<bin))
		{
			flag = p2_f(N,Q,dirname);
		}
		else
		{
			flag = p_f(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
	{
		if (Q==(1<<bin))
		{
			flag = p2_c(N,Q,dirname);
		}
		else
		{
			flag = p_c(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
	{
		if (Q==(1<<bin))
		{
			flag = p2_f_p(N,Q,dirname);
		}
		else
		{
			flag = p_f_p(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
	{
		if (Q==(1<<bin))
		{
			flag = p2_c_0p(N,Q,dirname);
		}
		else
		{
			flag = p_c_0p(N,Q,dirname);
		}
	}
}
else if (strcmp(option1,"pf")==0)
{
	if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
	{
		if (Q==(1<<bin))
		{
			flag = pf2_f(N,Q,dirname);
		}
		else
		{
			flag = pf_f(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
	{
		if (Q==(1<<bin))
		{
			flag = pf2_c(N,Q,dirname);
		}
		else
		{
			flag = pf_c(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
	{
		if (Q==(1<<bin))
		{
			flag = pf2_f_p(N,Q,dirname);
		}
		else
		{
			flag = pf_f_p(N,Q,dirname);
		}
	}
	else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
	{
		if (Q==(1<<bin))
		{
			flag = pf2_c_0p(N,Q,dirname);
		}
		else
		{
			flag = pf_c_0p(N,Q,dirname);
		}
	}
}

if (flag==0)
{
	time(&toc);
	totaltime = difftime(toc,tic);
	printf("\nThe total time in seconds was %gs.",totaltime);
}

return 0;
}




/*****************************************************/
/****************Full transfer matrices***************/
/*****************************************************/

/*******************************/
unsigned char i_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char xh;
char filename[256];

sprintf(filename,"%s/i_f_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	xh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d\n",(bit_sum(n^m)+xh));
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_c(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char xh;
char filename[256];

sprintf(filename,"%s/i_c_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	xh = bit_sum(n^circ_single_lshift(n,N));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d\n",(bit_sum(n^m)+xh));
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_f(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char xh,uh;
char filename[256];

sprintf(filename,"%s/if_f_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Free row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	xh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
	uh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+xh),(uh));
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_c(const unsigned char N, char* dirname)
{
unsigned long long n;
unsigned long long m;
FILE *fid;
unsigned char xh,uh;
char filename[256];

sprintf(filename,"%s/if_c_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}

/*Cylindrical row boundary conditions*/
for(n=0; n<(1ULL<<N); n++)
{
	xh = bit_sum(n^circ_single_lshift(n,N));
	uh = N-(bit_sum(n));
	for(m=0; m<(1ULL<<N); m++)
	{
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+xh),(uh));
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*****************************************************/
/**************Reduced transfer matrices**************/
/*****************************************************/

/*******************************/
unsigned char i_f_p(const unsigned char N, char* dirname)
{
unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char xmax=2*N-1;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((xmax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_f_p_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_f(N,&bitarray,&reflec,&total);
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
		xh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
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
						melement[(bit_sum(nn^mm)+xh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<xmax+1;p++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char i_c_0p(const unsigned char N, char* dirname)
{
unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char xmax=2*N;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char*) malloc((xmax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/i_c_0p_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_c(N,&bitarray,&reflec,&order,&total);
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
		xh = bit_sum(n^circ_single_lshift(n,N));
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
								melement[(bit_sum(nn^mm)+xh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<xmax+1;p++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_f_p(const unsigned char N, char* dirname)
{
unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,nn,mm;
unsigned char xmax=2*N-1;
unsigned char umax=N;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((xmax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N-1); n++)
{
	melement[n] = (unsigned char*) calloc(umax+1, sizeof(unsigned char));
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


sprintf(filename,"%s/if_f_p_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_f(N,&bitarray,&reflec,&total);
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
		xh = bit_sum((~1ULL&n)^(~1ULL&circ_single_lshift(n,N)));
		uh = N-(bit_sum(n));
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
						melement[(bit_sum(nn^mm)+xh)][(uh)]++;
						mm=bit_reflection(mm,N);
					}
					nn=bit_reflection(nn,N);
				}
				for (p=0;p<xmax+1;p++)
				{
					for (q=0;q<umax+1;q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
for (n=0ULL;n<(3*N-1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char if_c_0p(const unsigned char N, char* dirname)
{
unsigned char flag;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,nn,mm;
unsigned char xmax=2*N;
unsigned char umax=N;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,uh,rtotal=0ULL,ctotal=0ULL;

melement = (unsigned char**) malloc((xmax+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N+1); n++)
{
	melement[n] = (unsigned char*) calloc(umax+1, sizeof(unsigned char));
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

sprintf(filename,"%s/if_c_0p_%d.txt",dirname,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_c(N,&bitarray,&reflec,&order,&total);
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
		xh = bit_sum(n^circ_single_lshift(n,N));
		uh = N-(bit_sum(n));
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
								melement[((bit_sum(nn^mm))+xh)][(uh)]++;
								mm = circ_single_lshift(mm,N);
							}
							nn = circ_single_lshift(nn,N);
						}
						mm = bit_reflection(mm,N);
					}
					nn = bit_reflection(nn,N);
				}
				for (p=0;p<xmax+1;p++)
				{
					for (q=0;q<umax+1;q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
for (n=0ULL;n<(3*N+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}



/*****************************************************/
/********************Potts Models*********************/
/*****************************************************/

/*******************************/
unsigned char p2_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long xtest,xh,xinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	xtest = n ^ circ_bin_lshift(n,N,bin);
	xh=0ULL;
	for (p=1;p<N;p++)
	{
		xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}

	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{	
		xtest = n ^ m;
		xinter=0ULL;
		for (p=0;p<N;p++)
		{
			xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	
		fprintf(fid,"%llu\n",xh+xinter);
	}
}
	
printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p2_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long xtest,xh,xinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_c_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	xtest = n ^ circ_bin_lshift(n,N,bin);
	xh=0ULL;
	for (p=0;p<N;p++)
	{
		xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		xtest = n ^ m;
		xinter=0ULL;
		for (p=0;p<N;p++)
		{
			xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%llu\n",xh+xinter);
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long xtest,xh,xinter,uh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_f_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	xtest = n ^ circ_bin_lshift(n,N,bin);
	xh=0ULL;
	for (p=1;p<N;p++)
	{
		xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	uh = 0ULL;
	for (p=0;p<N;p++)
	{
		uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		xtest = n ^ m;
		xinter=0ULL;
		for (p=0;p<N;p++)
		{
			xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%llu %llu\n",xh+xinter,uh);
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long xtest,xh,xinter,uh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_c_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}


//Find transfer matrix
for (n=0ULL;n<(1ULL<<(bin*N));n++)
{
	xtest = n ^ circ_bin_lshift(n,N,bin);
	xh=0ULL;
	for (p=0;p<N;p++)
	{
		xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	uh = 0ULL;
	for (p=0;p<N;p++)
	{
		uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
	}
	
	for (m=0ULL;m<(1ULL<<(bin*N));m++)
	{
		xtest = n ^ m;
		xinter=0ULL;
		for (p=0;p<N;p++)
		{
			xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		fprintf(fid,"%llu %llu\n",xh+xinter,uh);
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long xtest,xh,xinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_f_%llu_%d.txt",dirname,Q,N);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
	
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				xtest = n ^ m;
				xinter=0ULL;
				for (p=0;p<N;p++)
				{
					xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
			
				fprintf(fid,"%llu\n",xh+xinter);
			}
		}
	}
}
	

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long xtest,xh,xinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_c_%llu_%d.txt",dirname,Q,N);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				xtest = n ^ m;
				xinter=0ULL;
				for (p=0;p<N;p++)
				{
					xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%llu\n",xh+xinter);
			}
		}
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long xtest,xh,xinter,uh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_f_%llu_%d.txt",dirname,Q,N);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				xtest = n ^ m;
				xinter=0ULL;
				for (p=0;p<N;p++)
				{
					xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%llu %llu\n",xh+xinter,uh);
			}
		}
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long xtest,xh,xinter,uh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_c_%llu_%d.txt",dirname,Q,N);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		
		for (m=0ULL;m<(1ULL<<(bin*N));m++)
		{
			bitfrac2=lldiv(m,csize);
			if (((pnums[bitfrac2.quot]&(1<<bitfrac2.rem))>>bitfrac2.rem)==0)
			{
				xtest = n ^ m;
				xinter=0ULL;
				for (p=0;p<N;p++)
				{
					xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				
				fprintf(fid,"%llu %llu\n",xh+xinter,uh);
			}
		}
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}



/****Reduced Potts transfer matrices****/

/*******************************/
unsigned char p2_f_p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}


melement = (unsigned char*) malloc((3*N-1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}


sprintf(filename,"%s/p_f_p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_bin_f(N,bin,&bitarray,&reflec,&total);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
						xtest = nn ^ mm;
						xinter=0ULL;
						for (r=0;r<N;r++)
						{
							xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[xh+xinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N-1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}

/*******************************/
unsigned char p2_c_0p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}

melement = (unsigned char*) malloc((3*N+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/p_c_0p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_bin_c(N,bin,&bitarray,&reflec,&order,&total);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
								xtest = nn ^ mm;
								xinter=0ULL;
								for (t=0;t<N;t++)
								{
									xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[xh+xinter]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N+1;p++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;

}

/*******************************/
unsigned char pf2_f_p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,uh,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}


melement = (unsigned char**) malloc((3*N-1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N-1); n++)
{
	melement[n] = (unsigned char*) calloc(2*N+1, sizeof(unsigned char));
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


sprintf(filename,"%s/pf_f_p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_bin_f(N,bin,&bitarray,&reflec,&total);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
						xtest = nn ^ mm;
						xinter=0ULL;
						for (r=0;r<N;r++)
						{
							xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[xh+xinter][uh]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N-1;p++)
				{
					for (q=0;q<(2*N+1);q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
for (n=0ULL;n<(3*N-1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char pf2_c_0p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,uh,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}

melement = (unsigned char**) malloc((3*N+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N+1); n++)
{
	melement[n] = (unsigned char*) calloc(2*N+1, sizeof(unsigned char));
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

sprintf(filename,"%s/pf_c_0p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = simple_red_bin_c(N,bin,&bitarray,&reflec,&order,&total);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
								xtest = nn ^ mm;
								xinter=0ULL;
								for (t=0;t<N;t++)
								{
									xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[xh+xinter][uh]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N+1;p++)
				{
					for (q=0;q<(2*N+1);q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
for (n=0ULL;n<(3*N+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char p_f_p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,nn,mm,sum;
unsigned char* pnums;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((3*N-1)*sizeof(unsigned char));
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


sprintf(filename,"%s/p_f_p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = arb_red_bin_f(N,bin,&bitarray,&reflec,&total,pnums);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
						xtest = nn ^ mm;
						xinter=0ULL;
						for (r=0;r<N;r++)
						{
							xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[xh+xinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N-1;p++)
				{
					if (melement[p]!=0)
					{
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
						melement[p]=0; //reset  - use memset once at end of loop? No, because this is more efficient, not all values are non-zero
					}
				}
			}
		}
		ctotal=0ULL;
	}
}

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char p_c_0p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((3*N+1)*sizeof(unsigned char));
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

sprintf(filename,"%s/p_c_0p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = arb_red_bin_c(N,bin,&bitarray,&reflec,&order,&total,pnums);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
								xtest = nn ^ mm;
								xinter=0ULL;
								for (t=0;t<N;t++)
								{
									xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[xh+xinter]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N+1;p++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char pf_f_p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,nn,mm,sum;
unsigned char* pnums;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,uh,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char**) malloc((3*N-1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N-1); n++)
{
	melement[n] = (unsigned char*) calloc(2*N+1, sizeof(unsigned char));
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


sprintf(filename,"%s/pf_f_p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = arb_red_bin_f(N,bin,&bitarray,&reflec,&total,pnums);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=1;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
						xtest = nn ^ mm;
						xinter=0ULL;
						for (r=0;r<N;r++)
						{
							xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*r))>>bin*r)==0ULL);
						}
						melement[xh+xinter][uh]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N-1;p++)
				{
					for (q=0;q<(2*N+1);q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
for (n=0ULL;n<(3*N-1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*******************************/
unsigned char pf_c_0p(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char flag;
unsigned char bin;
unsigned long long total=0ULL;
unsigned char* bitarray;
unsigned char* reflec;
unsigned char* order;
FILE* fid;
char filename[256];	
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long xh,xtest,xinter,uh,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char**) malloc((3*N+1)*sizeof(unsigned char*));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for(n = 0ULL; n < (3*N+1); n++)
{
	melement[n] = (unsigned char*) calloc(2*N+1, sizeof(unsigned char));
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

sprintf(filename,"%s/pf_c_0p_%llu_%d.txt",dirname,Q,N);
fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file. %s\n",strerror(errno));
	return 1;
}	

//Compute unique configurations
flag = arb_red_bin_c(N,bin,&bitarray,&reflec,&order,&total,pnums);
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
		xtest = n ^ circ_bin_lshift(n,N,bin);
		xh=0ULL;
		for (p=0;p<N;p++)
		{
			xh = xh + (((xtest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		uh = 0ULL;
		for (p=0;p<N;p++)
		{
			uh = uh + (((n & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
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
								xtest = nn ^ mm;
								xinter=0ULL;
								for (t=0;t<N;t++)
								{
									xinter = xinter + (((xtest & (((1ULL<<bin)-1ULL)<<bin*t))>>bin*t)==0ULL);
								}
								melement[xh+xinter][uh]++;
								mm = circ_bin_lshift(mm,N,bin);
							}
							nn = circ_bin_lshift(nn,N,bin);
						}
						mm = bit_reflection_bin(mm,N,bin);
					}
					nn = bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<3*N+1;p++)
				{
					for (q=0;q<(2*N+1);q++)
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

printf("\nFile  ../%s  created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)order);
for (n=0ULL;n<(3*N+1);n++)
{
	free((void*)melement[n]);
}
free((void*)melement);
return 0;
}


/*****************************************************/
/********************Row Reductions*******************/
/*****************************************************/

/*******************************/
unsigned char simple_red_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
{
unsigned long long n,m;
const unsigned char csize=8*sizeof(unsigned char);
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
unsigned char simple_red_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=8*sizeof(unsigned char);

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
unsigned char simple_red_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
{
unsigned long long n,m;
const unsigned char csize=8*sizeof(unsigned char);
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
unsigned char simple_red_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=8*sizeof(unsigned char);

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
unsigned char arb_red_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m;
const unsigned char csize=8*sizeof(unsigned char);
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
unsigned char arb_red_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=8*sizeof(unsigned char);

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
unsigned char arb_red_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m;
const unsigned char csize=8*sizeof(unsigned char);
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
unsigned char arb_red_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
{
unsigned long long n,m,p;
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
const unsigned char csize=8*sizeof(unsigned char);

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












