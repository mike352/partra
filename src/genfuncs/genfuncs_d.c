#include "partra_genfuncs.h"
#include "partra_genfuncs_d.h"

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
unsigned char ordering[2]={0,0}, remaining;
unsigned char flag;
double z1=0,z2=0;

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
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				//printf("n=%llu m=%llu m[%llu]=%hhu m[%llu]=%f\n",n,m,omsize[1]*p,imatrix[n][m][1][imsize[1]*p+ordering[1]],omsize[1]*p+1,pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*imatrix[n][m][1][imsize[1]*p+2]);
				(*omatrix)[n][m][1][omsize[1]*p]=imatrix[n][m][1][imsize[1]*p+ordering[1]];
				(*omatrix)[n][m][1][omsize[1]*p+1] = pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*imatrix[n][m][1][imsize[1]*p+2];
			}
		}
	}
}
else if ((remaining==0)&(imsize[1]==3)) //double substitution
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
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*pow(z2,imatrix[n][m][1][imsize[1]*p+ordering[1]])*imatrix[n][m][1][imsize[1]*p+2];
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
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
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
		(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
		for (p=0;p<imatrix[n][m][0][0];p++)
		{
			(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*imatrix[n][m][1][imsize[1]*p+1];
		}
	}
}

return 0;
}

/*******************************/
unsigned char matrix_fprintf_d(double**** matrix, unsigned long long* msize, char* filename, unsigned char prec)
{
unsigned long long n,m,p,q;
int fcheck=0;
char precformat[256];
FILE* fid;

sprintf(precformat," %%.%hhuf",prec);

fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file %s. %s\n",filename,strerror(errno));
	return 1;
}

for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		for (p=0ULL;p<(int)matrix[n][m][0][0];p++)
		{
			fcheck = fprintf(fid,"%llu %llu",n+1,m+1);
			if (fcheck<0)
			{
				printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
				fclose(fid);
				return 1;
			}
			for (q=0;q<msize[1]-1;q++)
			{
				fcheck = fprintf(fid," %hhu",(unsigned char) matrix[n][m][1][p*msize[1]+q]);
				if (fcheck<0)
				{
					printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
					fclose(fid);
					return 1;
				}
			}
			fcheck = fprintf(fid,precformat,matrix[n][m][1][p*msize[1]+q]);
			if (fcheck<0)
			{
				printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
				fclose(fid);
				return 1;
			}
			fprintf(fid,"\n");
		}
	}
}
printf("File %s created.\n",filename);
fclose(fid);
return 0;
}




/*************************************************/
/*************Double Complex Functions************/
/*************************************************/

/*******************************/
unsigned char matrix_alloc_dc(double complex***** matrix,const unsigned long long* msize, const unsigned char N)
{
unsigned long long n,m,p,q,r,s;	

*matrix = (double complex****) malloc(msize[0]*sizeof(double complex***)); 
if (*matrix==NULL)
{
	printf("\nERROR: Could not allocate memory.");
	return 2;
}
for (n=0ULL;n<msize[0];n++)
{
	(*matrix)[n] = (double complex***) malloc(msize[0]*sizeof(double complex**));
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
		(*matrix)[n][m] = (double complex**) malloc(2*sizeof(double complex*));
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
			(*matrix)[n][m][p] = (double complex*) calloc((msize[1]*N-2)*p+2,sizeof(double complex));
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
void matrix_free_dc(double complex**** matrix, const unsigned long long* msize)
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
unsigned char matrix_sub_dc(double complex***** omatrix, unsigned long long* omsize, unsigned char**** imatrix, unsigned long long* imsize, char* which, ...)
{
unsigned long long n,m,p;
unsigned char ordering[2]={0,0}, remaining;
unsigned char flag;
double complex z1=0,z2=0;

va_list vl;
va_start(vl,which);

if (strcmp(which,"u")==0)
{
	if (imsize[1]==3)
	{
		ordering[0]=0;
		ordering[1]=1;
		remaining=1;
		z1=va_arg(vl,double complex);
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		z1=va_arg(vl,double complex);
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
		z1=va_arg(vl,double complex);
		omsize[1]=imsize[1]-1;
	}
	else if (imsize[1]==2)
	{
		ordering[0]=0;
		remaining=0;
		z1=va_arg(vl,double complex);
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
	z1=va_arg(vl,double complex);
	z2=va_arg(vl,double complex);
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
	z1=va_arg(vl,double complex);
	z2=va_arg(vl,double complex);
	omsize[1]=imsize[1]-2;
}
else
{
	printf("ERROR: Incorrect choice of variable substitution. Only \"u\", \"x\", \"ux\", \"xu\" allowed.\n");
	return 3;
}

omsize[0]=imsize[0];
flag = matrix_alloc_dc(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
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
			if (cabs((*omatrix)[n][m][0][1])<cabs(imatrix[n][m][0][1]))
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double complex*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double complex));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_dc((*omatrix),omsize);
					return 2;
				}
			}
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				//printf("n=%llu m=%llu m[%llu]=%hhu m[%llu]=%f\n",n,m,omsize[1]*p,imatrix[n][m][1][imsize[1]*p+ordering[1]],omsize[1]*p+1,pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*imatrix[n][m][1][imsize[1]*p+2]);
				(*omatrix)[n][m][1][omsize[1]*p]=imatrix[n][m][1][imsize[1]*p+ordering[1]];
				(*omatrix)[n][m][1][omsize[1]*p+1] = pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*imatrix[n][m][1][imsize[1]*p+2];
			}
		}
	}
}
else if ((remaining==0)&(imsize[1]==3)) //double substitution
{
	for (n=0ULL;n<omsize[0];n++)
	{
		for (m=0ULL;m<omsize[0];m++)
		{
			if (cabs((*omatrix)[n][m][0][1])<cabs(imatrix[n][m][0][1]))
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double complex*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double complex));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_dc((*omatrix),omsize);
					return 2;
				}
			}
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
			for (p=0;p<imatrix[n][m][0][0];p++)
			{
				(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p+ordering[0]])*pow(z2,imatrix[n][m][1][imsize[1]*p+ordering[1]])*imatrix[n][m][1][imsize[1]*p+2];
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
			if (cabs((*omatrix)[n][m][0][1])<cabs(imatrix[n][m][0][1]))
			{
				(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
				(*omatrix)[n][m][1] = (double complex*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double complex));
				if ((*omatrix)[n][m][1]==NULL)
				{
					printf("ERROR: Unable to reallocate memory.\n");
					matrix_free_dc((*omatrix),omsize);
					return 2;
				}
			}
			(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
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
unsigned char matrix_sub_dc_dc(double complex***** omatrix, unsigned long long* omsize, double complex**** imatrix, unsigned long long* imsize, double complex z1)
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
flag = matrix_alloc_dc(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
if (flag!=0)
{
	return flag;
}

for (n=0ULL;n<omsize[0];n++)
{
	for (m=0ULL;m<omsize[0];m++)
	{
		if (cabs((*omatrix)[n][m][0][1])<cabs(imatrix[n][m][0][1]))
		{
			(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
			(*omatrix)[n][m][1] = (double complex*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double complex));
			if ((*omatrix)[n][m][1]==NULL)
			{
				printf("ERROR: Unable to reallocate memory.\n");
				matrix_free_dc((*omatrix),omsize);
				return 2;
			}
		}
		(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
		for (p=0;p<cabs(imatrix[n][m][0][0]);p++)
		{
			(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*imatrix[n][m][1][imsize[1]*p+1];
		}
	}
}

return 0;
}


/*******************************/
unsigned char matrix_sub_dc_d(double complex***** omatrix, unsigned long long* omsize, double**** imatrix, unsigned long long* imsize, double complex z1)
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
flag = matrix_alloc_dc(omatrix,omsize,imatrix[0][0][0][1]);  //allocate each row to previous row allocation
if (flag!=0)
{
	return flag;
}

for (n=0ULL;n<omsize[0];n++)
{
	for (m=0ULL;m<omsize[0];m++)
	{
		if (cabs((*omatrix)[n][m][0][1])<imatrix[n][m][0][1])
		{
			(*omatrix)[n][m][0][1] = imatrix[n][m][0][0];
			(*omatrix)[n][m][1] = (double complex*) realloc((*omatrix)[n][m][1],omsize[1]*imatrix[n][m][0][0]*sizeof(double complex));
			if ((*omatrix)[n][m][1]==NULL)
			{
				printf("ERROR: Unable to reallocate memory.\n");
				matrix_free_dc((*omatrix),omsize);
				return 2;
			}
		}
		(*omatrix)[n][m][0][0]=imatrix[n][m][0][0];
		for (p=0;p<imatrix[n][m][0][0];p++)
		{
			(*omatrix)[n][m][1][p] = pow(z1,imatrix[n][m][1][imsize[1]*p])*imatrix[n][m][1][imsize[1]*p+1];
		}
	}
}

return 0;
}


/*******************************/
unsigned char matrix_fprintf_dc(double complex**** matrix, unsigned long long* msize, char* filename, unsigned char prec)
{
unsigned long long n,m,p,q;
int fcheck=0;
char precformat[256];
FILE* fid;

sprintf(precformat," %%.%hhuf %%.%hhuf",prec,prec);

fid = fopen(filename,"w");
if (fid == NULL)
{
	printf("\nERROR: Could not create output file %s. %s\n",filename,strerror(errno));
	return 1;
}

for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		for (p=0ULL;p<(int)matrix[n][m][0][0];p++)
		{
			fcheck = fprintf(fid,"%llu %llu",n+1,m+1);
			if (fcheck<0)
			{
				printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
				fclose(fid);
				return 1;
			}
			for (q=0;q<msize[1]-1;q++)
			{
				fcheck = fprintf(fid," %hhu",(unsigned char) matrix[n][m][1][p*msize[1]+q]); //type-cast back to char
				if (fcheck<0)
				{
					printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
					fclose(fid);
					return 1;
				}
			}
			fcheck = fprintf(fid,precformat,creal(matrix[n][m][1][p*msize[1]+q]),cimag(matrix[n][m][1][p*msize[1]+q]));
			if (fcheck<0)
			{
				printf("ERROR: Problem writing to output file %s. %s\n",filename,strerror(errno));
				fclose(fid);
				return 1;
			}
			
			fprintf(fid,"\n");
		}
	}
}
printf("File %s created.\n",filename);
fclose(fid);
return 0;
}



