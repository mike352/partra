#include "partra_lapack.h"


int compeigs(const void* x, const void* y)
{
  complex double xx = *(complex double*)x;
  complex double yy = *(complex double*)y;
  if (cabs(xx) > cabs(yy)) return -1;
  if (cabs(xx) < cabs(yy)) return  1;
  return 0;
}


unsigned char eigs_lapack(double**** matrix, unsigned long long* msize)
{
unsigned long long n,m;

//Setting up LAPACK
int eignum=msize[0];
char JOBVL ='N';   // Compute Right eigenvectors
char JOBVR ='V';   // Do not compute Left eigenvectors
double complex VL[1];
int LDVL = 1; 
int LDVR = msize[0];
int LWORK = 4*msize[0]; 
int RRWORK = 2*msize[0];
double complex* WORK;
//double* WORK;
double complex* RWORK;
int INFO;
WORK = (double complex*)malloc(LWORK*sizeof(double complex));
//WORK = (double*)malloc(LWORK*sizeof(double));
if (WORK==NULL)
{
	printf("ERROR: Could not allocate memory.\n");
	return 2;
}
RWORK = (double complex*)malloc(RRWORK*sizeof(double complex));
if (RWORK==NULL)
{
	printf("ERROR: Could not allocate memory.\n");
	free(WORK);
	return 2;
}

int (*fp)(const void*,const void*) = compeigs;

double complex* M;
double complex* eigenvectors;
double complex* eigenvalues;

M = (double complex*) malloc((msize[0]*msize[0])*sizeof(double complex));
if (M==NULL)
{
	printf("ERROR: Could not allocate memory.\n");
	free(WORK);
	free(RWORK);
	return 2;
}
eigenvectors = (double complex*)malloc((msize[0]*msize[0]) * sizeof(double complex));
if (eigenvectors==NULL)
{
	printf("ERROR: Could not allocate memory.\n");
	free(WORK);
	free(RWORK);
	free(M);
	return 2;
}
eigenvalues = (double complex*)malloc((msize[0]) * sizeof(double complex));
if (eigenvalues==NULL)
{
	printf("ERROR: Could not allocate memory.\n");
	free(WORK);
	free(RWORK);
	free(M);
	free(eigenvectors);
	return 2;
}

for(n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		M[n*msize[0]+m] = matrix[n][m][1][0];
	}
}

//Use LAPACK to calculate all eigenvalues
zgeev_(&JOBVL, &JOBVR, &eignum, M, &eignum, eigenvalues, VL, &LDVL, eigenvectors, &LDVR, WORK, &LWORK, RWORK, &INFO);
//Sort the eigenvalues (they're not always sorted correctly!)
qsort(eigenvalues,msize[0],sizeof(eigenvalues[0]),fp);	

if (eignum<msize[0])
{
	printf("The number of eigenvalues %d found is less than the size of the matrix %llu.\n",eignum,msize[0]);
	return 5;
}

for (n=0;n<msize[0];n++)
{
	printf("%f %f\n",creal(eigenvalues[n]),cimag(eigenvalues[n]));
}

free(WORK);
free(RWORK);
free(M);
free(eigenvectors);
free(eigenvalues);
return 0;
}


