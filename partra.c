/* Options:
1. Ising or Ising in a field or Potts or Potts in a field
2. row boundary conditions 
3. full or reduced transfer matrix
4. row size
5. square or triangular lattice

Output file:
All data is output to a file in a (constructed) directory "data" from the current directory. If the OS is not unix, apple, or windows, it will ask for the directory to place the file. Filename is displayed on screen after successful completion of calculations and writing of data to file. 

Output file format legend: 
r = row
c = column
n: u^n = exp(-2E/kT)^n, energy exponent
m: x^m = exp(-2H/kT)^m, field exponent. In q-Potts, H-field only acts on q=1
A: A x^n u^m, constant 

Output file format:
1. Full transfer matrix:
n m
(r and c are assumed)
When H=0, m is not output

2. Reduced transfer matrix
r c A n m
When H=0, m is not output

Design principle:
Only unsigned char or unsigned long long data types are used for consistency. The sole exception currently is the dcheck variable which checks whether a directory was successfully created, which is of type int because output from mkdir commands can be -1. 

To Add:
1. Add api versions
2. Add field and/or energy substitutor 
3. Incorporate partra library into old equimodular curve finder program
4. Write a portable bit architecture tester 
5. Write a test suite


*/

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <complex.h>

//Only used to create a data directory or check that one already exists
 #if (defined(__unix__) || defined(__APPLE__))
	#include <sys/stat.h> 
	#define OS 1 
	int _mkdir(char*); //prototyping to make warning go away when compiling with gcc -Wall
#elif (defined(_WIN16) || defined(_WIN32)|| defined(_WIN64))
	#include <direct.h> 
	#define OS 0
#else //Otherwise ask
	#define OS -1 
#endif

typedef unsigned char**** Matrix;
typedef unsigned long long**** Matrix_ll;
typedef long double complex**** Matrix_ldc;
typedef unsigned char*** Row;
typedef unsigned long long*** Row_ll;
typedef long double complex*** Row_ldc;

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


//Transfer matrix functions
unsigned char i_sq_f(const unsigned char, char*); //Ising full transfer matrix, free row b.c.
unsigned char i_sq_c(const unsigned char, char*); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_sq_f(const unsigned char, char*); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_sq_c(const unsigned char, char*); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_sq_f_r(const unsigned char, char*); //Ising reduced transfer matrix, free row b.c.
unsigned char i_sq_c_r(const unsigned char, char*); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_sq_f_r(const unsigned char, char*); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_sq_c_r(const unsigned char, char*); //Ising in a field reduced transfer matrix, cylindrical row b.c.

unsigned char i_tri_f(const unsigned char, char*); //Ising full transfer matrix, free row b.c.
unsigned char i_tri_c(const unsigned char, char*); //Ising full transfer matrix, cylindrical row b.c.
unsigned char if_tri_f(const unsigned char, char*); //Ising in a field full transfer matrix, free row b.c.
unsigned char if_tri_c(const unsigned char, char*); //Ising in a field full transfer matrix, cylindrical row b.c.
unsigned char i_tri_f_r(const unsigned char, char*); //Ising reduced transfer matrix, free row b.c.
unsigned char i_tri_c_r(const unsigned char, char*); //Ising reduced transfer matrix, cylindrical row b.c.
unsigned char if_tri_f_r(const unsigned char, char*); //Ising in a field reduced transfer matrix, free row b.c.
unsigned char if_tri_c_r(const unsigned char, char*); //Ising in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_sq_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c.
unsigned char p_sq_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_sq_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p_tri_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c.
unsigned char p_tri_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c.
unsigned char pf_tri_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c.
unsigned char p_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c.
unsigned char p_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c.
unsigned char pf_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c.
unsigned char pf_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c.

unsigned char p2_sq_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_sq_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_sq_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

unsigned char p2_tri_f(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c(const unsigned char, const unsigned long long, char*); //Potts full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c(const unsigned char, const unsigned long long, char*); //Potts in a field full transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char p2_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char p2_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts reduced transfer matrix, cylindrical row b.c., when q is a power of 2
unsigned char pf2_tri_f_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, free row b.c., when q is a power of 2
unsigned char pf2_tri_c_r(const unsigned char, const unsigned long long, char*); //Potts in a field reduced transfer matrix, cylindrical row b.c., when q is a power of 2

//Reduction functions
unsigned char red_simple_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N
unsigned char red_simple_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N
unsigned char red_simple_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*); //cylindrical b.c. row reduction on all integers up to N*bin
unsigned char red_simple_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*); //free b.c. row reduction on all integers up to N*bin

unsigned char red_gen_c(const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_f(const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_bin_c(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //cylindrical b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array
unsigned char red_gen_bin_f(const unsigned char,const unsigned char,unsigned char**,unsigned char**,unsigned long long*,unsigned char*); //free b.c. row reduction on integers up to N*bin whose corresponding bits are 0 in supplied bit array


//General functions
unsigned char matrix_pow_ll(unsigned long long*****, unsigned long long*****, const unsigned long long*,const unsigned long long);

unsigned char matrix_alloc(unsigned char*****,const unsigned long long*,const unsigned char);
void matrix_free(unsigned char****,const unsigned long long*);
unsigned char row_alloc(unsigned char****,const unsigned long long*,const unsigned char);
void row_free(unsigned char***,const unsigned long long*);
unsigned char matrix_setadd(unsigned char*****, const unsigned long long*, const unsigned long long, const unsigned long long, const unsigned char*);
unsigned char row_setadd(unsigned char****, const unsigned long long*, const unsigned long long, const unsigned char*);

unsigned char matrix_alloc_ll(unsigned long long*****,const unsigned long long*,const unsigned char);
void matrix_free_ll(unsigned long long****,const unsigned long long*);
unsigned char row_alloc_ll(unsigned long long****,const unsigned long long*,const unsigned char);
void row_free_ll(unsigned long long***,const unsigned long long*);
unsigned char matrix_setadd_ll(unsigned long long*****, const unsigned long long*, const unsigned long long, const unsigned long long, const unsigned long long*);
unsigned char row_setadd_ll(unsigned long long****, const unsigned long long*, const unsigned long long, const unsigned long long*);

/*Size of Ising 0+ sector follows OEIS series A000029*/
/*Size of Ising 0 sector follows OEIS series A000031*/
/*Size of Ising + sector follows OEIS series A005418*/

int main(void)
{
unsigned char numoptions = 6;

unsigned char row_max_size=0;
FILE* fid=NULL;
unsigned char fcheck;
char* options[numoptions+1];
unsigned long long test;
unsigned char N=0;
char Qs[256];
unsigned long long Q=0ULL;
unsigned char bin=0;
char dirname[256];
int dcheck; //Needs to be of type int for use with mkdir
time_t tic;
time_t toc;
double totaltime;
char option0[256];
char option1[256];
char option2[256];
char option3[256];
char option4[256];
char option5[256];
char tmp[256];
unsigned char flag=~0;
unsigned char n,m;
unsigned char repeat=1;


//Discover max row size
test = 1ULL;
while (test!=0ULL)
{
	row_max_size++;
	test = test << 1;
}

//Create directory.
sprintf(dirname,"data");
if (OS) //UNIX
{
	dcheck = mkdir(dirname,S_IRWXU | S_IRWXG | S_IRWXO);
	if ((dcheck==-1)&(errno!=EEXIST))
	{
		printf("\nERROR: Could not create output directory. %s\n",strerror(errno));
		return 0;
	}
}
else if (OS==0) //Windows
{
	dcheck = _mkdir(dirname);
	if ((dcheck==-1)&(errno!=EEXIST))
	{
		printf("\nERROR: Could not create output directory. %s\n",strerror(errno));
		return 0;
	}
}
else if (OS==-1) //Ask instead
{
	printf("\nDirectory to place output files: ");
	scanf("%s",dirname);
}


//Ask for input file
printf("\nUse input file \"partra_setup.txt\"? (y,n): ");
scanf("%s",option0);
if ((strcmp(option0,"y")!=0)&(strcmp(option0,"n")))
{
	printf("\nERROR: Wrong input.");
	return 0;
}
if (strcmp(option0,"y")==0)
{
	fid = fopen("partra_setup.txt","r");
	if (fid==NULL)
	{
		printf("\nERROR: Could not open input file. %s\n",strerror(errno));
		return 0;
	}
}

//Repeating Loop
while ((strcmp(option0,"y")==0)&(repeat==1))
{
	if (strcmp(option0,"y")==0)
	{
		for (n=0;n<numoptions+1;n++)
		{
			options[n] = (char*) malloc(256*sizeof(char));
			if (options[n]==NULL)
			{
				printf("ERROR: Could not allocate memory.\n");
				for (m=0;m<n;m++)
				{
					free(options[m]);
				}
				fclose(fid);
				return 0;
			}
		}
		for (n=0;n<numoptions+1;n++)
		{
			fcheck = fscanf(fid,"%s",options[n]);
			if (fcheck!=1)
			{
				printf("ERROR: Problem reading option %d in the input file.",n+1);
				for (n=0;n<numoptions+1;n++)
				{
					free(options[n]);
				}
				fclose(fid);
				return 0;
			}
		}
		
		//Choice of Model
		strcpy(option1,options[0]);
		if ((strcmp(option1,"i")!=0)&(strcmp(option1,"if")!=0)&(strcmp(option1,"p")!=0)&(strcmp(option1,"pf")!=0))
		{
			printf("\nERROR: Wrong model input value.");
			flag=1;
		}
		
		//Value of q (for Potts only)
		if ((strcmp(option1,"p")==0)|(strcmp(option1,"pf")==0))
		{
			if (strspn(options[1],"1234567890")<strlen(options[1]))
			{
				printf("\nERROR: Invalid q value."); //prevents negative values which would get reinterpreted, floats, and others
				flag=1;
			}
			else
			{
				Q=strtoll(options[1],NULL,10);
				if ((Q==LONG_MAX)|(Q==LONG_MIN))
				{
					printf("\nERROR: Given q value is too large. %s.",strerror(errno));
					flag=1;
				}
				else if (Q==0)
				{
					printf("\nERROR: Invalid q value."); //problem converting string Qs to unsigned long long Q
					flag=1;
				}
				else if (Q<3ULL)
				{
					printf("\nERROR: q should be greater than 2.");
					flag=1;
				}
				else if (Q>=3ULL)
				{
					while((1ULL<<bin)<Q)
					{
						bin++; //calculate bit array bin size
					}
				}
			}
		}
		
		//Choice of row boundary condition
		strcpy(option2,options[2]);
		if ((strcmp(option2,"f")!=0)&(strcmp(option2,"c")!=0))
		{
			printf("\nERROR: Wrong row boundary condition input.");
			flag=1;
		}
		
		//Choice of full or reduced matrix
		strcpy(option3,options[3]);
		if ((strcmp(option3,"f")!=0)&(strcmp(option3,"r")!=0))
		{
			printf("\nERROR: Wrong input for choice of full or reduced matrix.");
			flag=1;
		}
		
		//Size of row
		if (strspn(options[4],"1234567890")<strlen(options[4]))
		{
			printf("\nERROR: Invalid N value"); //prevents negative values which would get reinterpreted, floats, and others
			flag=1;
		}
		else
		{
			N=strtoll(options[4],NULL,10);
			if ((strcmp(option1,"i")==0)|(strcmp(option1,"if")==0))
			{
				if (N<1)
				{
					printf("\nERROR: Row size should be greater than 0.\n");
					flag=1;
				}
				else if (N>row_max_size)
				{
					printf("\nERROR: Your machine can only do %d-bit calcuations.\n       Limit row size to %d.\n",row_max_size,row_max_size);
					flag=1;
				}
			}
			else if ((strcmp(option1,"p")==0)|(strcmp(option1,"pf")==0))
			{
				if (N<1)
				{
					printf("\nERROR: Row size should be greater than 0.\n");
					flag=1;
				}
				else if (N>row_max_size/bin)
				{
					printf("\nERROR: Your machine can only do %d-bit calcuations.\n       Limit row size to %d.\n",row_max_size,row_max_size/(bin*N));
					flag=1;
				}
			}
		}
		
		//Choice of lattice
		strcpy(option4,options[5]);
		if ((strcmp(option4,"s")!=0)&(strcmp(option4,"t")!=0))
		{
			printf("\nERROR: Wrong lattice type input.");
			flag=1;
		}
		
		//Repeat program
		strcpy(option5,options[6]);
		if ((strcmp(option5,"y")!=0)&(strcmp(option5,"n")!=0))
		{
			printf("\nERROR: Wrong input of whether to continue reading input file.");
			flag=1;
		}
		if (strcmp(option5,"n")==0) //last loop
		{
			repeat=0;
			for (n=0;n<numoptions+1;n++)
			{
				free(options[n]);
				
			}
			fclose(fid);
		}
		else 
		{
			if (fgets(tmp,2,fid)==NULL) //Skip a blank line
			{
				printf("\nERROR: Problem reading input file after option %d.",numoptions+1);
				flag=1;
			}			
		}
		
		//Exit if any of the options were bad
		if (flag==1)
		{
			for (n=0;n<numoptions+1;n++)
			{
				free(options[n]);
				
			}
			fclose(fid);
			return 0;
		}
		
	}
	else
	{
		//Model choice
		printf("Ising or Ising in a field or Potts or Potts in a field  (i,if,p,pf): ");
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

		//Lattice
		printf("Square or triangular lattice (s,t): ");
		scanf("%s",option4);
		if ((strcmp(option4,"s")!=0)&(strcmp(option4,"t")!=0))
		{
			printf("\nERROR: Wrong input.");
			return 0;
		}
	}

	//Choose a function
	time(&tic);
	if (strcmp(option4,"s")==0)
	{
		if (strcmp(option1,"i")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				flag = i_sq_f(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				flag = i_sq_c(N,dirname);
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				flag = i_sq_f_r(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				flag = i_sq_c_r(N,dirname);
			}
		}
		else if (strcmp(option1,"if")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				flag = if_sq_f(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				flag = if_sq_c(N,dirname);
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				flag = if_sq_f_r(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				flag = if_sq_c_r(N,dirname);
			}
		}
		else if (strcmp(option1,"p")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_sq_f(N,Q,dirname);
				}
				else
				{
					flag = p_sq_f(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_sq_c(N,Q,dirname);
				}
				else
				{
					flag = p_sq_c(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_sq_f_r(N,Q,dirname);
				}
				else
				{
					flag = p_sq_f_r(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_sq_c_r(N,Q,dirname);
				}
				else
				{
					flag = p_sq_c_r(N,Q,dirname);
				}
			}
		}
		else if (strcmp(option1,"pf")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_sq_f(N,Q,dirname);
				}
				else
				{
					flag = pf_sq_f(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_sq_c(N,Q,dirname);
				}
				else
				{
					flag = pf_sq_c(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_sq_f_r(N,Q,dirname);
				}
				else
				{
					flag = pf_sq_f_r(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_sq_c_r(N,Q,dirname);
				}
				else
				{
					flag = pf_sq_c_r(N,Q,dirname);
				}
			}
		}
	}
	else if (strcmp(option4,"t")==0)
	{
		if (strcmp(option1,"i")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				flag = i_tri_f(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				flag = i_tri_c(N,dirname);
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				flag = i_tri_f_r(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				flag = i_tri_c_r(N,dirname);
			}
		}
		else if (strcmp(option1,"if")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				flag = if_tri_f(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				flag = if_tri_c(N,dirname);
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				flag = if_tri_f_r(N,dirname);
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				flag = if_tri_c_r(N,dirname);
			}
		}
		else if (strcmp(option1,"p")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_tri_f(N,Q,dirname);
				}
				else
				{
					flag = p_tri_f(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_tri_c(N,Q,dirname);
				}
				else
				{
					flag = p_tri_c(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_tri_f_r(N,Q,dirname);
				}
				else
				{
					flag = p_tri_f_r(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = p2_tri_c_r(N,Q,dirname);
				}
				else
				{
					flag = p_tri_c_r(N,Q,dirname);
				}
			}
		}
		else if (strcmp(option1,"pf")==0)
		{
			if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_tri_f(N,Q,dirname);
				}
				else
				{
					flag = pf_tri_f(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_tri_c(N,Q,dirname);
				}
				else
				{
					flag = pf_tri_c(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_tri_f_r(N,Q,dirname);
				}
				else
				{
					flag = pf_tri_f_r(N,Q,dirname);
				}
			}
			else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
			{
				if (Q==(1<<bin))
				{
					flag = pf2_tri_c_r(N,Q,dirname);
				}
				else
				{
					flag = pf_tri_c_r(N,Q,dirname);
				}
			}
		}
	}
if (flag==0)
{
	time(&toc);
	totaltime = difftime(toc,tic);
	printf("\n   The total time in seconds was %gs.",totaltime);
}
}

return 0;
}




/*****************************************************/
/**********Full Ising square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char i_sq_f(const unsigned char N, char* dirname)
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
unsigned char i_sq_c(const unsigned char N, char* dirname)
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
unsigned char if_sq_f(const unsigned char N, char* dirname)
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
unsigned char if_sq_c(const unsigned char N, char* dirname)
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
unsigned char i_tri_f(const unsigned char N, char* dirname)
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
	    fprintf(fid,"%d\n",(bit_sum(n^m)+bit_sum(n^(~1ULL&circ_single_lshift(m,N)))+uh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char i_tri_c(const unsigned char N, char* dirname)
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
unsigned char if_tri_f(const unsigned char N, char* dirname)
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
	    fprintf(fid,"%d %d\n",(bit_sum(n^m)+bit_sum(n^(~1ULL&circ_single_lshift(m,N)))+uh),(xh));
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char if_tri_c(const unsigned char N, char* dirname)
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
						melement[(bit_sum(nn^mm)+bit_sum(nn^(~1ULL&circ_single_lshift(mm,N)))+uh)]++;
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
						melement[(bit_sum(nn^mm)+bit_sum(nn^(~1ULL&circ_single_lshift(mm,N)))+uh)][(xh)]++;
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



/*****************************************************/
/*********Full Potts square transfer matrices*********/
/*****************************************************/

/*******************************/
unsigned char p2_sq_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_%llu_%d.txt",dirname,Q,N);
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
	
		fprintf(fid,"%llu\n",uh+uinter);
	}
}
	
printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p2_sq_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_%llu_%d.txt",dirname,Q,N);
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
		
		fprintf(fid,"%llu\n",uh+uinter);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_%llu_%d.txt",dirname,Q,N);
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
		
		fprintf(fid,"%llu %llu\n",uh+uinter,xh);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_sq_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_%llu_%d.txt",dirname,Q,N);
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
		
		fprintf(fid,"%llu %llu\n",uh+uinter,xh);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p_sq_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_f_%llu_%d.txt",dirname,Q,N);
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
			
				fprintf(fid,"%llu\n",uh+uinter);
			}
		}
	}
}
	

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_sq_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_sq_c_%llu_%d.txt",dirname,Q,N);
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
				
				fprintf(fid,"%llu\n",uh+uinter);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_f_%llu_%d.txt",dirname,Q,N);
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
				
				fprintf(fid,"%llu %llu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_sq_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_sq_c_%llu_%d.txt",dirname,Q,N);
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
				
				fprintf(fid,"%llu %llu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}



/*****************************************************/
/*******Full Potts triangular transfer matrices*******/
/*****************************************************/

/*******************************/
unsigned char p2_tri_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_tri_f_%llu_%d.txt",dirname,Q,N);
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
	for (p=1;p<N;p++) //p starts at 1 for free b.c.
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
		utest = n ^ circ_bin_lshift(m,N,bin);
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		fprintf(fid,"%llu\n",uh+uinter);
	}
}
	
printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p2_tri_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_tri_c_%llu_%d.txt",dirname,Q,N);
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
		utest = n ^ circ_bin_lshift(m,N,bin); 
		for (p=0;p<N;p++) 
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		fprintf(fid,"%llu\n",uh+uinter);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_tri_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_tri_f_%llu_%d.txt",dirname,Q,N);
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
	for (p=1;p<N;p++) //p starts at 1 for free b.c.
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
		utest = n ^ circ_bin_lshift(m,N,bin); 
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		fprintf(fid,"%llu %llu\n",uh+uinter,xh);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char pf2_tri_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_tri_c_%llu_%d.txt",dirname,Q,N);
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
		utest = n ^ circ_bin_lshift(m,N,bin); 
		for (p=0;p<N;p++)
		{
			uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
		}
		fprintf(fid,"%llu %llu\n",uh+uinter,xh);
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
return 0;
}


/*******************************/
unsigned char p_tri_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_tri_f_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
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
				utest = n ^ circ_bin_lshift(m,N,bin);
				for (p=1;p<N;p++) //p starts at 1 for free b.c.
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				fprintf(fid,"%llu\n",uh+uinter);
			}
		}
	}
}
	

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char p_tri_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/p_tri_c_%llu_%d.txt",dirname,Q,N);
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
				utest = n ^ circ_bin_lshift(m,N,bin);
				for (p=0;p<N;p++) 
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				fprintf(fid,"%llu\n",uh+uinter);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_tri_f(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_tri_f_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++) //p starts at 1 for free b.c.
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
				utest = n ^  circ_bin_lshift(m,N,bin);
				for (p=1;p<N;p++) //p starts at 1 for free b.c.
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				fprintf(fid,"%llu %llu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}


/*******************************/
unsigned char pf_tri_c(const unsigned char N, const unsigned long long Q, char* dirname)
{
const unsigned char csize=8*sizeof(unsigned char);
FILE* fid;
char filename[256];	
unsigned char bin=0;
unsigned long long n,m,p,sum;
unsigned char* pnums;
lldiv_t bitfrac,bitfrac2;
unsigned long long utest,uh,uinter,xh;

while((1ULL<<bin)<Q)
{
	bin++;
}

sprintf(filename,"%s/pf_tri_c_%llu_%d.txt",dirname,Q,N);
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
				utest = n ^ circ_bin_lshift(m,N,bin);
				for (p=0;p<N;p++)
				{
					uinter = uinter + (((utest & (((1ULL<<bin)-1ULL)<<bin*p))>>bin*p)==0ULL);
				}
				fprintf(fid,"%llu %llu\n",uh+uinter,xh);
			}
		}
	}
}

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)pnums);
return 0;
}



/*****************************************************/
/*******Reduced Potts square transfer matrices********/
/*****************************************************/

/*******************************/
unsigned char p2_sq_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N-1;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}


melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}


sprintf(filename,"%s/p_sq_f_r_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++)
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
						melement[uh+uinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}

/*******************************/
unsigned char p2_sq_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
if ((melement==NULL))
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

sprintf(filename,"%s/p_sq_c_r_%llu_%d.txt",dirname,Q,N);
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
unsigned char pf2_sq_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N-1;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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


sprintf(filename,"%s/pf_sq_f_r_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++)
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
unsigned char pf2_sq_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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

sprintf(filename,"%s/pf_sq_c_r_%llu_%d.txt",dirname,Q,N);
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


/*******************************/
unsigned char p_sq_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N-1;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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


sprintf(filename,"%s/p_sq_f_r_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++)
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
						melement[uh+uinter]++;
						mm=bit_reflection_bin(mm,N,bin);
					}
					nn=bit_reflection_bin(nn,N,bin);
				}
				for (p=0;p<umax+1;p++)
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

printf("\nFile  ../%s created.",filename);
fclose(fid);
free((void*)bitarray);
free((void*)reflec);
free((void*)melement);
return 0;
}


/*******************************/
unsigned char p_sq_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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

sprintf(filename,"%s/p_sq_c_r_%llu_%d.txt",dirname,Q,N);
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
unsigned char pf_sq_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N-1;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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


sprintf(filename,"%s/pf_sq_f_r_%llu_%d.txt",dirname,Q,N);
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
		for (p=1;p<N;p++)
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
unsigned char pf_sq_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 2*N;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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

sprintf(filename,"%s/pf_sq_c_r_%llu_%d.txt",dirname,Q,N);
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
/*****Reduced Potts triangular transfer matrices******/
/*****************************************************/

/*******************************/
unsigned char p2_tri_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 3*N-2;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}


melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
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
unsigned char p2_tri_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
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
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

while((1ULL<<bin)<Q)
{
	bin++;
}

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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
unsigned char pf2_tri_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 3*N-2;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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
unsigned char pf2_tri_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
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
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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


/*******************************/
unsigned char p_tri_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 3*N-2;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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
						fprintf(fid,"%llu %llu %llu %llu\n",rtotal,ctotal,melement[p]/(1ULL+flip2),p); //add normalization constant only to row vector
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
unsigned char p_tri_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
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
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char* melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,utest,uinter,rtotal=0ULL,ctotal=0ULL;

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

melement = (unsigned char*) malloc((umax+1)*sizeof(unsigned char));
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
unsigned char pf_tri_f_r(const unsigned char N, const unsigned long long Q, char* dirname)
{
unsigned char umax = 3*N-2;
unsigned char xmax = N;

unsigned char flag;
unsigned char bin=0;
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
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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
unsigned char pf_tri_c_r(const unsigned char N, const unsigned long long Q, char* dirname)
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
const unsigned char csize=8*sizeof(unsigned char);
unsigned long long n,m,p,q,r,s,t,nn,mm,sum;
unsigned char* pnums;
unsigned char** melement; 
lldiv_t bitfrac, bitfrac2, bitfrac3, bitfrac4;
unsigned char flip, flip2;
unsigned long long uh,utest,uinter,xh,rtotal=0ULL,ctotal=0ULL;

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
/********************Row Reductions*******************/
/*****************************************************/

/*******************************/
unsigned char red_simple_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
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
unsigned char red_simple_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
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
unsigned char red_simple_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total)
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
unsigned char red_simple_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total)
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
unsigned char red_gen_f(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
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
unsigned char red_gen_c(const unsigned char N,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
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
unsigned char red_gen_bin_f(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned long long* total,unsigned char *numbers)
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
unsigned char red_gen_bin_c(const unsigned char N,const unsigned char bin,unsigned char **bitarray, unsigned char **reflec, unsigned char** order, unsigned long long* total,unsigned char *numbers)
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






/*****************************************************/
/*******************General functions*****************/
/*****************************************************/



/*******************************/
unsigned char matrix_alloc(unsigned char***** matrix,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,m,p,q,r,s;	

*matrix = (unsigned char****) malloc(msize[0]*sizeof(unsigned char***)); 
if (*matrix==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*matrix)[n] = (unsigned char***) malloc(msize[0]*sizeof(unsigned char**));
	if (((*matrix)[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (q=0ULL;q<n;q++)
		{
			free((void*)(*matrix)[q]);
		}
		return 2;
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		(*matrix)[n][m] = (unsigned char**) malloc(2*sizeof(unsigned char*));
		if (((*matrix)[n][m]==NULL))
		{
			printf("\nERROR: Could not allocate memory.");
			for (q=0ULL;q<n;q++)
			{
				for (r=0ULL;r<m;r++)
				{
					free((void*)(*matrix)[q][r]);
				}
				free((void*)(*matrix)[q]);
			}
			return 2;
		}
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		for (p=0ULL;p<2;p++)
		{
			(*matrix)[n][m][p] = (unsigned char*) calloc((msize[1]*N-2)*p+2,sizeof(unsigned char));
			if (((*matrix)[n][m][p]==NULL))
			{
				printf("\nERROR: Could not allocate memory.");
				for (q=0ULL;q<n;q++)
				{
					for (r=0ULL;r<m;r++)
					{
						for (s=0ULL;s<p;s++)
						{
							free((void*)(*matrix)[q][r][s]);
						}
						free((void*)(*matrix)[q][r]);
					}
					free((void*)(*matrix)[q]);
				}
				return 2;
			}
		}
		(*matrix)[n][m][0][1]=N; //initial size of matrix[n][m][1][] is msize[1]*N
	}
}


return 0;
}


/*******************************/
void matrix_free(unsigned char**** matrix, const unsigned long long* msize)
{
unsigned long long q,r,s;

for (q=0ULL;q<msize[0];q++)
{
	for (r=0ULL;r<msize[0];r++)
	{
		for (s=0ULL;s<2;s++)
		{
			free((void*)matrix[q][r][s]);
		}
		free((void*)matrix[q][r]);
	}
	free((void*)matrix[q]);
}
free((void*)matrix);

}


/*******************************/
unsigned char row_alloc(unsigned char**** row,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,p,q,s;	

*row = (unsigned char***) malloc(msize[0]*sizeof(unsigned char**)); 
if (*row==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*row)[n] = (unsigned char**) malloc(2*sizeof(unsigned char*));
	if (((*row)[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (q=0ULL;q<n;q++)
		{
			free((void*)(*row)[q]);
		}
		return 2;
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (p=0ULL;p<2;p++)
	{
		(*row)[n][p] = (unsigned char*) calloc((msize[1]*N-2)*p+2,sizeof(unsigned char));
		if (((*row)[n][p]==NULL))
		{
			printf("\nERROR: Could not allocate memory.");
			for (q=0ULL;q<n;q++)
			{
				for (s=0ULL;s<p;s++)
				{
					free((void*)(*row)[q][s]);
				}
				free((void*)(*row)[q]);
			}
			return 2;
		}
	}
	(*row)[n][0][1]=N; //initial size of row[n][1][] is msize[1]*N
}


return 0;
}


/*******************************/
void row_free(unsigned char*** row, const unsigned long long* msize)
{
unsigned long long q,s;

for (q=0ULL;q<msize[0];q++)
{
	for (s=0ULL;s<2;s++)
	{
		free((void*)row[q][s]);
	}
	free((void*)row[q]);
}
free((void*)row);

}


/*******************************/
unsigned char matrix_setadd(unsigned char***** matrix, const unsigned long long* msize, const unsigned long long n, const unsigned long long m, const unsigned char* v)
{
unsigned long long q,r;
unsigned char test=1;


//Search and set appropriately
for (q=0ULL;q<(*matrix)[n][m][0][0];q++)
{
	for (r=0ULL;r<msize[1]-1;r++)
	{
		test = test & ((*matrix)[n][m][1][msize[1]*q+r]==v[r]);
	}
	if (test)
	{
		(*matrix)[n][m][1][msize[1]*q+msize[1]-1] = (*matrix)[n][m][1][msize[1]*q+msize[1]-1] + v[msize[1]-1];
		return 0;
	}
}

//Couldn't find it:

//There's still room allocated
if ((*matrix)[n][m][0][0]<(*matrix)[n][m][0][1])
{
	for (r=0ULL;r<msize[1];r++)
	{
		(*matrix)[n][m][1][msize[1]*(*matrix)[n][m][0][0]+r]=v[r];
	}
	(*matrix)[n][m][0][0]++;
}
else //Not enough room - re-allocate
{		
	//(*matrix)[n][m][1] = (unsigned char*) realloc((*matrix)[n][m][1],(msize[1]*(*matrix)[n][m][0][1]+msize[1]*msize[0])*sizeof(unsigned char)); //add another total
	(*matrix)[n][m][1] = (unsigned char*) realloc((*matrix)[n][m][1],(msize[1]*(*matrix)[n][m][0][1]+msize[1]*(*matrix)[n][m][0][1]*(*matrix)[n][m][0][1])*sizeof(unsigned char));	//add square of current amount
	if ((*matrix)[n][m][1] != NULL)
	{
		//(*matrix)[n][m][0][1] = (*matrix)[n][m][0][1]+msize[0]; //added another total
		(*matrix)[n][m][0][1] = (*matrix)[n][m][0][1]+(*matrix)[n][m][0][1]*(*matrix)[n][m][0][1]; //added square of current amount
		for (r=0ULL;r<msize[1];r++)
		{
			(*matrix)[n][m][1][msize[1]*(*matrix)[n][m][0][0]+r]=v[r];
		}
		(*matrix)[n][m][0][0]++;
	}
	else
	{
		printf("\nERROR: Could not re-allocate memory.");
		return 2;
	}
}

return 0;
}


/*******************************/
unsigned char row_setadd(unsigned char**** row, const unsigned long long* msize, const unsigned long long n, const unsigned char* v)
{
unsigned long long q,r;
unsigned char test=1;


//Search and set appropriately
for (q=0ULL;q<(*row)[n][0][0];q++)
{
	for (r=0ULL;r<msize[1]-1;r++)
	{
		test = test & ((*row)[n][1][msize[1]*q+r]==v[r]);
	}
	if (test)
	{
		(*row)[n][1][msize[1]*q+msize[1]-1] = (*row)[n][1][msize[1]*q+msize[1]-1] + v[msize[1]-1];
		return 0;
	}
}


//Couldn't find it:

//There's still room allocated
if ((*row)[n][0][0]<(*row)[n][0][1])
{
	for (r=0ULL;r<msize[1];r++)
	{
		(*row)[n][1][msize[1]*(*row)[n][0][0]+r]=v[r];
	}
	(*row)[n][0][0]++;
}
else //Not enough room - re-allocate
{		
	//(*row)[n][1] = (unsigned char*) realloc((*row)[n][1],(msize[1]*(*row)[n][0][1]+msize[1]*msize[0])*sizeof(unsigned char)); //add another total
	(*row)[n][1] = (unsigned char*) realloc((*row)[n][1],(msize[1]*(*row)[n][0][1]+msize[1]*(*row)[n][0][1]*(*row)[n][0][1])*sizeof(unsigned char));	//add square of current amount
	if ((*row)[n][1] != NULL)
	{
		//(*row)[n][0][1] = (*row)[n][0][1]+msize[0]; //added another total
		(*row)[n][0][1] = (*row)[n][0][1]+(*row)[n][0][1]*(*row)[n][0][1]; //added square of current amount
		for (r=0ULL;r<msize[1];r++)
		{
			(*row)[n][1][msize[1]*(*row)[n][0][0]+r]=v[r];
		}
		(*row)[n][0][0]++;
	}
	else
	{
		printf("\nERROR: Could not re-allocate memory.");
		return 2;
	}
}

return 0;
}



//Note: this function is not efficient. 
/*******************************/
unsigned char matrix_pow_ll(unsigned long long***** fmatrix, unsigned long long***** matrix, const unsigned long long* msize,const unsigned long long M)
{
unsigned char flag;
unsigned long long n,q,r,s,t,u,v;
unsigned long long*** row;
unsigned long long* arrsend;

arrsend = malloc(msize[1]*sizeof(unsigned long long));
if (arrsend==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}

flag = row_alloc_ll(&row,msize,(*matrix)[0][0][0][1]);
if (flag!=0)
{
	free(arrsend);
	return flag;
}

flag = matrix_alloc_ll(fmatrix,msize,(*matrix)[0][0][0][1]);
if (flag!=0)
{
	free(arrsend);
	row_free_ll(row,msize);
	return flag;
}

//Copy matrix over to fmatrix
for (q=0ULL;q<msize[0];q++)
{
	for (r=0ULL;r<msize[0];r++)
	{
		(*fmatrix)[q][r][0][0]=(*matrix)[q][r][0][0];
		for (s=0ULL;s<msize[1]*(*matrix)[q][r][0][0];s++)
		{
			(*fmatrix)[q][r][1][s]=(*matrix)[q][r][1][s];
		}
	}
}


//If M=1, returns
//Otherwise, do multiplications
for (n=1ULL;n<M;n++)
{
	for (q=0ULL;q<msize[0];q++)
	{
		for (r=0ULL;r<msize[0];r++)
		{
			for (s=0ULL;s<msize[0];s++)
			{
				for (u=0ULL;u<(*fmatrix)[q][s][0][0];u++)
				{
					for (v=0ULL;v<(*matrix)[s][r][0][0];v++)
					{
						for (t=0ULL;t<msize[1]-1;t++)
						{
							arrsend[t]=(*fmatrix)[q][s][1][msize[1]*u+t]+(*matrix)[s][r][1][msize[1]*v+t];
						}
						arrsend[msize[1]-1] = ((*fmatrix)[q][s][1][msize[1]*u+msize[1]-1])*((*matrix)[s][r][1][msize[1]*v+msize[1]-1]);
						flag = row_setadd_ll(&row,msize,r,arrsend);
						if (flag!=0)
						{
							free(arrsend);
							matrix_free_ll(*matrix,msize);
							matrix_free_ll(*fmatrix,msize);
							row_free_ll(row,msize);
							return flag;
						}
					}
				}
			}
		}
		for (r=0ULL;r<msize[0];r++)
		{//printf("r=%d ",r);
			(*fmatrix)[q][r][0][0] = row[r][0][0];
			(*fmatrix)[q][r][0][1] = row[r][0][0]; //Make as big as number of valid entries
			(*fmatrix)[q][r][1] = (unsigned long long*) realloc((*fmatrix)[q][r][1],msize[1]*row[r][0][0]*sizeof(unsigned long long));
			
			//(*fmatrix)[q][r][0][1] = row[r][0][1]; //Make as big as row[r][1]
			//(*fmatrix)[q][r][1] = (unsigned long long*) realloc((*fmatrix)[q][r][1],msize[1]*row[r][0][1]*sizeof(unsigned long long));
			
			if ((*fmatrix)[q][r][1]!=NULL)
			{
				memcpy((*fmatrix)[q][r][1],row[r][1],msize[1]*row[r][0][0]*sizeof(unsigned long long)); //only copy valid entries
				
				//memcpy((*fmatrix)[q][r][1],row[r][1],msize[1]*row[r][0][1]*sizeof(unsigned long long)); //copy all of row[r][1]
				
				for (s=0ULL;s<msize[1]*row[r][0][0];s++)
				{
					row[r][1][s]=0ULL; //reset the non-zero values of row[r][1]
				}
				row[r][0][0]=0ULL;
			}
			else
			{
				printf("\nERROR: Could not allocate memory.");
				free(arrsend);
				matrix_free_ll(*matrix,msize);
				matrix_free_ll(*fmatrix,msize);
				row_free_ll(row,msize);
				return 2;
			}
		}
	}
}



row_free_ll(row,msize);
free(arrsend);
return 0;
}


/*******************************/
unsigned char matrix_alloc_ll(unsigned long long***** matrix,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,m,p,q,r,s;	

*matrix = (unsigned long long****) malloc(msize[0]*sizeof(unsigned long long***)); 
if (*matrix==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*matrix)[n] = (unsigned long long***) malloc(msize[0]*sizeof(unsigned long long**));
	if (((*matrix)[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (q=0ULL;q<n;q++)
		{
			free((void*)(*matrix)[q]);
		}
		return 2;
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		(*matrix)[n][m] = (unsigned long long**) malloc(2*sizeof(unsigned long long*));
		if (((*matrix)[n][m]==NULL))
		{
			printf("\nERROR: Could not allocate memory.");
			for (q=0ULL;q<n;q++)
			{
				for (r=0ULL;r<m;r++)
				{
					free((void*)(*matrix)[q][r]);
				}
				free((void*)(*matrix)[q]);
			}
			return 2;
		}
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		for (p=0ULL;p<2;p++)
		{
			(*matrix)[n][m][p] = (unsigned long long*) calloc((msize[1]*N-2)*p+2,sizeof(unsigned long long));
			if (((*matrix)[n][m][p]==NULL))
			{
				printf("\nERROR: Could not allocate memory.");
				for (q=0ULL;q<n;q++)
				{
					for (r=0ULL;r<m;r++)
					{
						for (s=0ULL;s<p;s++)
						{
							free((void*)(*matrix)[q][r][s]);
						}
						free((void*)(*matrix)[q][r]);
					}
					free((void*)(*matrix)[q]);
				}
				return 2;
			}
		}
		(*matrix)[n][m][0][1]=N; //initial size of matrix[n][m][1][] is msize[1]*N
	}
}


return 0;
}


/*******************************/
void matrix_free_ll(unsigned long long**** matrix, const unsigned long long* msize)
{
unsigned long long q,r,s;

for (q=0ULL;q<msize[0];q++)
{
	for (r=0ULL;r<msize[0];r++)
	{
		for (s=0ULL;s<2;s++)
		{
			free((void*)matrix[q][r][s]);
		}
		free((void*)matrix[q][r]);
	}
	free((void*)matrix[q]);
}
free((void*)matrix);

}


/*******************************/
unsigned char row_alloc_ll(unsigned long long**** row,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,p,q,s;	

*row = (unsigned long long***) malloc(msize[0]*sizeof(unsigned long long**)); 
if (*row==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*row)[n] = (unsigned long long**) malloc(2*sizeof(unsigned long long*));
	if (((*row)[n]==NULL))
	{
		printf("\nERROR: Could not allocate memory.");
		for (q=0ULL;q<n;q++)
		{
			free((void*)(*row)[q]);
		}
		return 2;
	}
}
for (n=0ULL;n<msize[0];n++)
{
	for (p=0ULL;p<2;p++)
	{
		(*row)[n][p] = (unsigned long long*) calloc((msize[1]*N-2)*p+2,sizeof(unsigned long long));
		if (((*row)[n][p]==NULL))
		{
			printf("\nERROR: Could not allocate memory.");
			for (q=0ULL;q<n;q++)
			{
				for (s=0ULL;s<p;s++)
				{
					free((void*)(*row)[q][s]);
				}
				free((void*)(*row)[q]);
			}
			return 2;
		}
	}
	(*row)[n][0][1]=N; //initial size of row[n][1][] is msize[1]*N
}


return 0;
}


/*******************************/
void row_free_ll(unsigned long long*** row, const unsigned long long* msize)
{
unsigned long long q,s;

for (q=0ULL;q<msize[0];q++)
{
	for (s=0ULL;s<2;s++)
	{
		free((void*)row[q][s]);
	}
	free((void*)row[q]);
}
free((void*)row);

}


/*******************************/
unsigned char matrix_setadd_ll(unsigned long long***** matrix, const unsigned long long* msize, const unsigned long long n, const unsigned long long m, const unsigned long long* v)
{
unsigned long long q,r;
unsigned char test=1;


//Search and set appropriately
for (q=0ULL;q<(*matrix)[n][m][0][0];q++)
{
	for (r=0ULL;r<msize[1]-1;r++)
	{
		test = test & ((*matrix)[n][m][1][msize[1]*q+r]==v[r]);
	}
	if (test)
	{
		(*matrix)[n][m][1][msize[1]*q+msize[1]-1] = (*matrix)[n][m][1][msize[1]*q+msize[1]-1] + v[msize[1]-1];
		return 0;
	}
}

//Couldn't find it:

//There's still room allocated
if ((*matrix)[n][m][0][0]<(*matrix)[n][m][0][1])
{
	for (r=0ULL;r<msize[1];r++)
	{
		(*matrix)[n][m][1][msize[1]*(*matrix)[n][m][0][0]+r]=v[r];
	}
	(*matrix)[n][m][0][0]++;
}
else //Not enough room - re-allocate
{		
	//(*matrix)[n][m][1] = (unsigned long long*) realloc((*matrix)[n][m][1],(msize[1]*(*matrix)[n][m][0][1]+msize[1]*msize[0])*sizeof(unsigned long long)); //add another total
	(*matrix)[n][m][1] = (unsigned long long*) realloc((*matrix)[n][m][1],(msize[1]*(*matrix)[n][m][0][1]+msize[1]*(*matrix)[n][m][0][1]*(*matrix)[n][m][0][1])*sizeof(unsigned long long));	//add square of current amount
	if ((*matrix)[n][m][1] != NULL)
	{
		//(*matrix)[n][m][0][1] = (*matrix)[n][m][0][1]+msize[0]; //added another total
		(*matrix)[n][m][0][1] = (*matrix)[n][m][0][1]+(*matrix)[n][m][0][1]*(*matrix)[n][m][0][1]; //added square of current amount
		for (r=0ULL;r<msize[1];r++)
		{
			(*matrix)[n][m][1][msize[1]*(*matrix)[n][m][0][0]+r]=v[r];
		}
		(*matrix)[n][m][0][0]++;
	}
	else
	{
		printf("\nERROR: Could not re-allocate memory.");
		return 2;
	}
}

return 0;
}


/*******************************/
unsigned char row_setadd_ll(unsigned long long**** row, const unsigned long long* msize, const unsigned long long n, const unsigned long long* v)
{
unsigned long long q,r;
unsigned char test=1;


//Search and set appropriately
for (q=0ULL;q<(*row)[n][0][0];q++)
{
	for (r=0ULL;r<msize[1]-1;r++)
	{
		test = test & ((*row)[n][1][msize[1]*q+r]==v[r]);
	}
	if (test)
	{
		(*row)[n][1][msize[1]*q+msize[1]-1] = (*row)[n][1][msize[1]*q+msize[1]-1] + v[msize[1]-1];
		return 0;
	}
}

//Couldn't find it:

//There's still room allocated
if ((*row)[n][0][0]<(*row)[n][0][1])
{
	for (r=0ULL;r<msize[1];r++)
	{
		(*row)[n][1][msize[1]*(*row)[n][0][0]+r]=v[r];
	}
	(*row)[n][0][0]++;
}
else //Not enough room - re-allocate
{		
	//(*row)[n][1] = (unsigned long long*) realloc((*row)[n][1],(msize[1]*(*row)[n][0][1]+msize[1]*msize[0])*sizeof(unsigned long long)); //add another total
	(*row)[n][1] = (unsigned long long*) realloc((*row)[n][1],(msize[1]*(*row)[n][0][1]+msize[1]*(*row)[n][0][1]*(*row)[n][0][1])*sizeof(unsigned long long));	//add square of current amount
	if ((*row)[n][1] != NULL)
	{
		//(*row)[n][0][1] = (*row)[n][0][1]+msize[0]; //added another total
		(*row)[n][0][1] = (*row)[n][0][1]+(*row)[n][0][1]*(*row)[n][0][1]; //added square of current amount
		for (r=0ULL;r<msize[1];r++)
		{
			(*row)[n][1][msize[1]*(*row)[n][0][0]+r]=v[r];
		}
		(*row)[n][0][0]++;
	}
	else
	{
		printf("\nERROR: Could not re-allocate memory.");
		return 2;
	}
}

return 0;
}

