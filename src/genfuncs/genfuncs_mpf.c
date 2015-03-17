#include "partra_genfuncs.h"
#include "partra_genfuncs_mpf.h"
#include "gmp.h"

//NOTE: Only integer exponents can be used. For non-integer exponents, MPFR will be needed to carry out powers. 

/*******************************/
unsigned char matrix_sub_d_mpf(double***** omatrix, unsigned long long* omsize, unsigned char**** imatrix, unsigned long long* imsize, char* which, ...)
{
unsigned long long n,m,p;
unsigned char ordering[2], remaining;
unsigned char flag;
mpf_t z1,z2,tmp1,tmp2;
unsigned char prec=21;

va_list vl;
va_start(vl,which);

if (strcmp(which,"u")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=0;
		ordering[1]=1;
		remaining=1;
		mpf_init2(z1,prec);
		mpf_init2(z2,prec);
		mpf_set_d(z1,va_arg(vl,double));
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		mpf_init2(z1,prec);
		mpf_init2(z2,prec);
		mpf_set_d(z1,va_arg(vl,double));
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
		mpf_init2(z1,prec);
		mpf_init2(z2,prec);
		mpf_set_d(z1,va_arg(vl,double));
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		mpf_init2(z1,prec);
		mpf_init2(z2,prec);
		mpf_set_d(z1,va_arg(vl,double));
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
	mpf_init2(z1,prec);
	mpf_init2(z2,prec);
	mpf_set_d(z1,va_arg(vl,double));
	mpf_set_d(z2,va_arg(vl,double));
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
	mpf_init2(z1,prec);
	mpf_init2(z2,prec);
	mpf_set_d(z2,va_arg(vl,double));
	mpf_set_d(z1,va_arg(vl,double));
	omsize[1]=imsize[1]-2;
}
else
{
	printf("ERROR: Incorrect choice of variable substitution. Only \"u\", \"x\", \"ux\", \"xu\" allowed.\n");
	return 3;
}
va_end(vl);

omsize[0]=imsize[0];
if (remaining==1)
{
	flag = matrix_alloc_d(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
}
else if (remaining==0)
{
	flag = matrix_alloc_d(omatrix,omsize,1);
}
if (flag!=0)
{
	mpf_clears(z1,z2,NULL);
	return flag;
}


mpf_init2(tmp1,prec);
mpf_init2(tmp2,prec);

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
					mpf_clears(z1,z2,tmp1,tmp2,NULL);
					matrix_free_d((*omatrix),omsize);
					return 2;
				}
			}
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][omsize[1]*p]=imatrix[n][m][1][imsize[1]*p+ordering[1]];
				mpf_pow_ui(tmp1,z1,imatrix[n][m][1][imsize[1]*p+ordering[0]]);
				mpf_mul_ui(tmp1,tmp1,imatrix[n][m][1][imsize[1]*p+2]);
				(*omatrix)[n][m][1][omsize[1]*p+1] = mpf_get_d(tmp1);
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
			(*omatrix)[n][m][0][0]=1;
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				mpf_pow_ui(tmp1,z1,imatrix[n][m][1][imsize[1]*p]);
				mpf_pow_ui(tmp2,z2,imatrix[n][m][1][imsize[1]*p+1]);
				mpf_mul_ui(tmp1,tmp1,imatrix[n][m][1][imsize[1]*p+2]);
				mpf_mul(tmp1,tmp1,tmp2);
				(*omatrix)[n][m][1][0] = (*omatrix)[n][m][1][0] + mpf_get_d(tmp1);
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
			(*omatrix)[n][m][0][0]=1;
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				mpf_pow_ui(tmp1,z1,imatrix[n][m][1][imsize[1]*p]);
				mpf_mul_ui(tmp1,tmp1,imatrix[n][m][1][imsize[1]*p+1]);
				(*omatrix)[n][m][1][0] = (*omatrix)[n][m][1][0] + mpf_get_d(tmp1);
			}
		}
	}
}

mpf_clears(z1,z2,tmp1,tmp2,NULL);
return 0;
}


/*******************************/
unsigned char matrix_sub_d_d_mpf(double***** omatrix, unsigned long long* omsize, double**** imatrix, unsigned long long* imsize, double t)
{
unsigned long long n,m,p;
unsigned char flag;
mpf_t z1,tmp1,tmp2;
unsigned char prec=21;

if (imsize[1]==2)
{
	mpf_init2(z1,prec);
	mpf_set_d(z1,t);
	omsize[1]=imsize[1]-1;
}
else
{
	printf("ERROR: No free variables to use for substitution\n");
	return 3;
}


omsize[0]=imsize[0];
flag = matrix_alloc_d(omatrix,omsize,1);  //allocate each row to previous row allocation
if (flag!=0)
{
	mpf_clear(z1);
	return flag;
}

mpf_init2(tmp1,prec);
mpf_init2(tmp2,prec);
for (n=0ULL;n<omsize[0];n++)
{
	for (m=0ULL;m<omsize[0];m++)
	{
		(*omatrix)[n][m][0][0]=1;
		for (p=0;p<imatrix[n][m][0][0];p++)
		{
			mpf_pow_ui(tmp1,z1,(unsigned long int)imatrix[n][m][1][imsize[1]*p]);
			mpf_set_d(tmp2,imatrix[n][m][1][imsize[1]*p+1]);
			mpf_mul(tmp1,tmp1,tmp2);
			(*omatrix)[n][m][1][0] = (*omatrix)[n][m][1][0] + mpf_get_d(tmp1);
		}
	}
}

mpf_clears(z1,tmp1,tmp2,NULL);
return 0;
}

