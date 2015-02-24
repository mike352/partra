#include "partra_genfuncs.h"
#include "gmp.h"

/*
There is a problem with this function. It needs to create a new matrix of a different type than the input matrix. 
*/

unsigned char matrix_sub_mpf(double***** omatrix, unsigned long long* omsize, unsigned char**** imatrix, unsigned long long* imsize, char* which, ...)
{
unsigned long long n,m,p,q;
unsigned char ordering[2], remaining;
mpf_t z,tmp;
mpf_init2(z,21);
mpf_init2(tmp,21);

va_list vl;
va_start(vl,which);
mpf_set_d(z,va_arg(vl,double));

if (strcmp(which,"u")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=0;
		ordering[1]=1;
		remaining=1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
	}
	else
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
	omsize[1]=imsize[1]-1;
}
else if (strcmp(which,"x")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=1;
		ordering[1]=0;
		remaining=1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
	}
	else
	{
		printf("ERROR: No free variables to use for substitution\n");
		return 3;
	}
	omsize[1]=imsize[1]-1;
}
else if (strcmp(which,"ux")==0)
{
	ordering[0]=0;
	ordering[1]=1;
	remaining=0;
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
	omsize[1]=imsize[1]-2;
}
else if (strcmp(which,"xu")==0)
{
	ordering[0]=1;
	ordering[1]=0;
	remaining=0;
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
	omsize[1]=imsize[1]-2;
}
else
{
	printf("ERROR: Incorrect choice of variable substitution. Only \"u\", \"x\", \"ux\", \"xu\" allowed.\n");
	return 3;
}

omsize[0]=imsize[0];
matrix_alloc_d(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation

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
				(*omatrix)[n][m][1][omsize[1]*p]=imatrix[n][m][1][imsize[1]*p+ordering[0]];
				mpf_pow_ui(tmp,z,imatrix[n][m][1][imsize[1]*p+ordering[1]]);
				mpf_mul_ui(tmp,tmp,imatrix[n][m][1][imsize[1]*p+2]);
				(*omatrix)[n][m][1][omsize[1]*p+1] = mpf_get_d(tmp);
			}
		}
	}
}
else if (remaining==0)
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
				for (q=0;q<imsize[1];q++) //continue working on this part 
				{
					mpf_pow_ui(tmp,z,imatrix[n][m][1][imsize[1]*p+ordering[q]]);
					mpf_mul_ui(tmp,tmp,imatrix[n][m][1][imsize[1]*p+2]);
				}
				(*omatrix)[n][m][1][p] = mpf_get_d(tmp);
			}
		}
	}
}

mpf_clear(tmp);
return 0;
}