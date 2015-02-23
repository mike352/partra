#include "partra_genfuncs.h"
#include "gmp.h"

/*
There is a problem with this function. It needs to create a new matrix of a different type than the input matrix. 
*/

unsigned char matrix_sub_mpf(unsigned char***** matrix, unsigned long long* msize, char* which, ...)
{
unsigned long long n,m,p;
unsigned char a,b;
mpf_t z,tmp;
mpf_init2(z,21);
mpf_init2(tmp,21);

va_list vl;
va_start(vl,which);
mpf_set_d(z,va_arg(vl,double));

if (strcmp(which,"u")==0)
{
	a=0;b=1;
}
else if (strcmp(which,"x")==0)
{
	a=1;b=0;
}
else
{
	printf("ERROR: Incorrect choice of variable substitution\n");
	return 3;
}


unsigned long long rowsize=(*matrix)[0][0][0][0];
double* row;
row = calloc(2*rowsize,sizeof(double));

for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		if ((*matrix)[n][m][0][0]>rowsize)
		{
			rowsize = (*matrix)[n][m][0][0];
			row = realloc(row,2*rowsize*sizeof(double));
			if (row==NULL)
			{
				printf("ERROR: Unable to reallocate memory.\n");
				free(row);
				return 2;
			}
		}
		for (p=0;p<(*matrix)[n][m][0][0];p++)
		{
			row[2*p]=(*matrix)[n][m][1][msize[1]*p+a];
			row[2*p+2] = mpf_get_d(mpf_mul_ui(tmp,mpf_pow_ui(tmp,z,(*matrix)[n][m][1][msize[1]*p+b]),(*matrix)[n][m][1][msize[1]*p+2])));
		}
		matrix[n][m][1] = realloc(matrix[n][m][1],2*(*matrix)[n][m][0][0]*sizeof(double));
		if (matrix[n][m][1] == NULL)
		{
			printf("ERROR: Unable to reallocate memory.\n");
			free(row);
			return 2;
		}
	}
}
msize[1]=2;

mpf_clear(tmp);
free(row);
return 0;
}