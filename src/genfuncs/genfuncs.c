#include "partra_genfuncs.h"
#include <math.h>

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


unsigned char matrix_fprintf(unsigned char***** matrix, unsigned long long* msize, char* filename, char* dirname)
{
unsigned long long n,m,p,q;
int fcheck=0;
char outfile[256];
FILE* fid;

sprintf(outfile,"%s/%s",dirname,filename);
fid = fopen(outfile,"w");
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		for (p=0;p<(*matrix)[n][m][0][0];p++)
		{
			fprintf(fid,"%llu %llu",n,m);
			for (q=0;q<msize[1];q++)
			{
				fcheck = fprintf(fid," %hhu",(*matrix)[n][m][1][p*msize[1]+q]);
				if (fcheck<0)
				{
					printf("ERROR: Problem writing to output file. %s\n",strerror(errno));
					return 1;
				}
			}
			fprintf(fid,"\n");
		}
	}
}
fclose(fid);



return 0;
}



/*******************************/
unsigned char matrix_alloc_d(double***** matrix,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,m,p,q,r,s;	

*matrix = (double****) malloc(msize[0]*sizeof(double***)); 
if (*matrix==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*matrix)[n] = (double***) malloc(msize[0]*sizeof(double**));
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
		(*matrix)[n][m] = (double**) malloc(2*sizeof(double*));
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
			(*matrix)[n][m][p] = (double*) calloc((msize[1]*N-2)*p+2,sizeof(double));
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
void matrix_free_d(double**** matrix, const unsigned long long* msize)
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
unsigned char matrix_sub_d(double***** omatrix, unsigned long long* omsize, unsigned char**** imatrix, unsigned long long* imsize, char* which, ...)
{
unsigned long long n,m,p;
unsigned char ordering[2], remaining;
unsigned char flag;
double z1,z2;

va_list vl;
va_start(vl,which);

if (strcmp(which,"u")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=0;
		ordering[1]=1;
		remaining=1;
		z1=va_arg(vl,double);
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		z1=va_arg(vl,double);
		omsize[1]=imsize[1]-1;
	}
	else
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
}
else if (strcmp(which,"x")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=1;
		ordering[1]=0;
		remaining=1;
		z1=va_arg(vl,double);
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		z1=va_arg(vl,double);
		omsize[1]=imsize[1]-1;
	}
	else
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
}
else if (strcmp(which,"ux")==0)
{
	if (imsize[1]==2)
	{
		printf("ERROR: Not enough free variables to use for substitution\n");
		return 3;
	}
	else if (imsize[1]<2)
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
	ordering[0]=0;
	ordering[1]=1;
	remaining=0;
	z1=va_arg(vl,double);
	z2=va_arg(vl,double);
	omsize[1]=imsize[1]-2;
}
else if (strcmp(which,"xu")==0)
{
	if (imsize[1]==2)
	{
		printf("ERROR: Not enough free variables to use for substitution\n");
		return 3;
	}
	else if (imsize[1]<2)
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
	ordering[0]=1;
	ordering[1]=0;
	remaining=0;
	z1=va_arg(vl,double);
	z2=va_arg(vl,double);
	omsize[1]=imsize[1]-2;
}
else
{
	printf("ERROR: Incorrect choice of variable substitution. Only \"u\", \"x\", \"ux\", \"xu\" allowed.\n");
	return 3;
}

omsize[0]=imsize[0];
flag = matrix_alloc_d(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
if (flag!=0)
{
	return flag;
}

if (remaining==1)
{
	for (n=0ULL;n<omsize[0];n++)
	{
		for (m=0ULL;m<omsize[0];m++)
		{
			if ((*omatrix)[n][m][0][1]<imatrix[n][m][0][1])
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_d((*omatrix),omsize);
					return 2;
				}
			}
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][omsize[1]*p]=imatrix[n][m][1][imsize[1]*p+ordering[1]];
				(*omatrix)[n][m][1][omsize[1]*p+1] = pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*imatrix[n][m][1][imsize[1]*p+2];
			}
		}
	}
}
else if ((remaining==0)&(imsize[1]==3))
{
	for (n=0ULL;n<omsize[0];n++)
	{
		for (m=0ULL;m<omsize[0];m++)
		{
			if ((*omatrix)[n][m][0][1]<imatrix[n][m][0][1])
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_d((*omatrix),omsize);
					return 2;
				}
			}
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*pow(z2,imatrix[n][m][1][imsize[1]*p+1])*imatrix[n][m][1][imsize[1]*p+2];
			}
		}
	}
}
else if ((remaining==0)&(imsize[1]==2))
{
	for (n=0ULL;n<omsize[0];n++)
	{
		for (m=0ULL;m<omsize[0];m++)
		{
			if ((*omatrix)[n][m][0][1]<imatrix[n][m][0][1])
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_d((*omatrix),omsize);
					return 2;
				}
			}
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*imatrix[n][m][1][imsize[1]*p+1];
			}
		}
	}
}

return 0;
}

/*******************************/
unsigned char matrix_sub_d_d(double***** omatrix, unsigned long long* omsize, double**** imatrix, unsigned long long* imsize, double z1)
{
unsigned long long n,m,p;
unsigned char flag;

if (imsize[1]==2)
{
	omsize[1]=imsize[1]-1;
}
else
{
	printf("ERROR: No free variables to use for substitution\n");
	return 3;
}


omsize[0]=imsize[0];
flag = matrix_alloc_d(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
if (flag!=0)
{
	return flag;
}

for (n=0ULL;n<omsize[0];n++)
{
	for (m=0ULL;m<omsize[0];m++)
	{
		if ((*omatrix)[n][m][0][1]<imatrix[n][m][0][1])
		{
			(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
			(*omatrix)[n][m][1] = (double*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double));
			if ((*omatrix)[n][m][1]==NULL)
			{
				printf("ERROR: Unable to reallocate memory.\n");
				matrix_free_d((*omatrix),omsize);
				return 2;
			}
		}
		for (p=0;p<imatrix[n][m][0][0];p++)
		{
			(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*imatrix[n][m][1][imsize[1]*p+1];
		}
	}
}

return 0;
}