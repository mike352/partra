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
Only unsigned char or unsigned long long data types are used for integers for consistency. The sole exception currently is the dcheck variable which checks whether a directory was successfully created, which is of type int because output from mkdir commands can be -1. For strings, only char arrays of size 256 are used for consistency.
*/

#include "partra.h"

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


//Selector function
unsigned char selector(const unsigned char, const unsigned long long, const unsigned char, const char*,const char*,const char*,const char*,const char*);

int main(int argc, char* argv[])
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
int dcheck,scheck; //Needs to be of type int for use with mkdir and for checking scanf
time_t tic;
time_t toc;
double totaltime;
char option1[256];
char option2[256];
char option3[256];
char option4[256];
char option5[256];
char tmp[256];
unsigned char flag=~0;
unsigned char n,m;
unsigned char repeat=0;


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
	scheck=scanf("%s",dirname);
	if (scheck<=0)
	{
		printf("ERROR: %s",strerror(errno));
		return 0;
	}
}


//Repeating loop of inputs
if (argc==2)
{
	fid = fopen(argv[1],"r");
	if (fid==NULL)
	{
		printf("\nERROR: Could not open input file. %s\n",strerror(errno));
		return 0;
	}
	repeat=1;
	while (repeat==1)
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
		//Choose a function
		time(&tic);
		flag = selector(N,Q,bin,dirname,option1,option2,option3,option4);
		if (flag==0)
		{
			time(&toc);
			totaltime = difftime(toc,tic);
			printf("\n   The total time in seconds was %gs.\n",totaltime);
		}
		else
		{
			return 0;
		}
	}
}
else if (argc==1) //Command line options
{
	//Model choice
	printf("Ising or Ising in a field or Potts or Potts in a field  (i,if,p,pf): ");
	scheck=scanf("%s",option1);
	if (scheck<=0)
	{
		printf("ERROR: %s",strerror(errno));
		return 0;
	}
	if ((strcmp(option1,"i")!=0)&(strcmp(option1,"if")!=0)&(strcmp(option1,"p")!=0)&(strcmp(option1,"pf")!=0))
	{
		printf("\nERROR: Wrong input.");
		return 0;
	}

	//Input q value, checking carefully, since very large numbers are allowed. 
	if ((strcmp(option1,"p")==0)|(strcmp(option1,"pf")==0))
	{
		printf("Potts q value  (3 to 2^%d-1): ",row_max_size);
		scheck=scanf("%255s",Qs);
		if (scheck<=0)
		{
			printf("ERROR: %s",strerror(errno));
			return 0;
		}
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
	scheck=scanf("%s",option2);
	if (scheck<=0)
	{
		printf("ERROR: %s",strerror(errno));
		return 0;
	}
	if ((strcmp(option2,"f")!=0)&(strcmp(option2,"c")!=0))
	{
		printf("\nERROR: Wrong input.");
		return 0;
	}

	//Ask for transfer matrix reduction
	printf("Full or reduced transfer matrix (f,r): ");
	scheck=scanf("%s",option3);
	if (scheck<=0)
	{
		printf("ERROR: %s",strerror(errno));
		return 0;
	}
	if ((strcmp(option3,"f")!=0)&(strcmp(option3,"r")!=0))
	{
		printf("\nERROR: Wrong input.");
		return 0;
	}

	//Ask for row size
	if ((strcmp(option1,"i")==0)|(strcmp(option1,"if")==0))
	{
		printf("Row size (1 to %hhu): ",row_max_size); 
		scheck=scanf("%hhu",&N);
		if (scheck<=0)
		{
			printf("ERROR: %s",strerror(errno));
			return 0;
		}
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
		scheck=scanf("%hhu",&N);
		if (scheck<=0)
		{
			printf("ERROR: %s",strerror(errno));
			return 0;
		}
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
	scheck=scanf("%s",option4);
	if (scheck<=0)
	{
		printf("ERROR: %s",strerror(errno));
		return 0;
	}
	if ((strcmp(option4,"s")!=0)&(strcmp(option4,"t")!=0))
	{
		printf("\nERROR: Wrong input.");
		return 0;
	}
	//Choose a function
	time(&tic);
	flag = selector(N,Q,bin,dirname,option1,option2,option3,option4);
	if (flag==0)
	{
		time(&toc);
		totaltime = difftime(toc,tic);
		printf("\n   The total time in seconds was %gs.\n",totaltime);
	}
}
else
{
	printf("\nERROR: Too many inputs.");
}
	

return 0;
}


// Selector function
unsigned char selector(const unsigned char N, const unsigned long long Q, const unsigned char bin, const char* dirname, const char* option1, const char* option2, const char* option3, const char* option4)
{
unsigned char flag=~0;

if (strcmp(option4,"s")==0)
{
	if (strcmp(option1,"i")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			flag = i_sq_f_f_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			flag = i_sq_c_f_file(N,dirname);
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			flag = i_sq_f_r_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			flag = i_sq_c_r_file(N,dirname);
		}
	}
	else if (strcmp(option1,"if")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			flag = if_sq_f_f_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			flag = if_sq_c_f_file(N,dirname);
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			flag = if_sq_f_r_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			flag = if_sq_c_r_file(N,dirname);
		}
	}
	else if (strcmp(option1,"p")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_sq_f_f_file(N,Q,dirname);
			}
			else
			{
				flag = p_sq_f_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_sq_c_f_file(N,Q,dirname);
			}
			else
			{
				flag = p_sq_c_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_sq_f_r_file(N,Q,dirname);
			}
			else
			{
				flag = p_sq_f_r_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_sq_c_r_file(N,Q,dirname);
			}
			else
			{
				flag = p_sq_c_r_file(N,Q,dirname);
			}
		}
	}
	else if (strcmp(option1,"pf")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_sq_f_f_file(N,Q,dirname);
			}
			else
			{
				flag = pf_sq_f_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_sq_c_f_file(N,Q,dirname);
			}
			else
			{
				flag = pf_sq_c_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_sq_f_r_file(N,Q,dirname);
			}
			else
			{
				flag = pf_sq_f_r_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_sq_c_r_file(N,Q,dirname);
			}
			else
			{
				flag = pf_sq_c_r_file(N,Q,dirname);
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
			flag = i_tri_f_f_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			flag = i_tri_c_f_file(N,dirname);
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			flag = i_tri_f_r_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			flag = i_tri_c_r_file(N,dirname);
		}
	}
	else if (strcmp(option1,"if")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			flag = if_tri_f_f_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			flag = if_tri_c_f_file(N,dirname);
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			flag = if_tri_f_r_file(N,dirname);
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			flag = if_tri_c_r_file(N,dirname);
		}
	}
	else if (strcmp(option1,"p")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_tri_f_f_file(N,Q,dirname);
			}
			else
			{
				flag = p_tri_f_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_tri_c_f_file(N,Q,dirname);
			}
			else
			{
				flag = p_tri_c_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_tri_f_r_file(N,Q,dirname);
			}
			else
			{
				flag = p_tri_f_r_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = p2_tri_c_r_file(N,Q,dirname);
			}
			else
			{
				flag = p_tri_c_r_file(N,Q,dirname);
			}
		}
	}
	else if (strcmp(option1,"pf")==0)
	{
		if ((strcmp(option2,"f")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_tri_f_f_file(N,Q,dirname);
			}
			else
			{
				flag = pf_tri_f_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"f")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_tri_c_f_file(N,Q,dirname);
			}
			else
			{
				flag = pf_tri_c_f_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"f")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_tri_f_r_file(N,Q,dirname);
			}
			else
			{
				flag = pf_tri_f_r_file(N,Q,dirname);
			}
		}
		else if ((strcmp(option2,"c")==0) & (strcmp(option3,"r")==0))
		{
			if (Q==(1<<bin))
			{
				flag = pf2_tri_c_r_file(N,Q,dirname);
			}
			else
			{
				flag = pf_tri_c_r_file(N,Q,dirname);
			}
		}
	}
}

return flag;
}