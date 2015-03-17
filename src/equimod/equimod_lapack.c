#include<stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES //In order to use M_PI for pi
#include<math.h>
#include<complex.h>


//......................................................................................................
//  Computes the eigenvalues of input matrix A, sorts them by magnitude, 
//  and computes the relative ratio of the top two
//.....................................................................................................
double EigenRatio( double complex *eigenvectorsVR, double complex *eigenvaluesW, double r,int *S, int *z, int N,int (*fp)(const void*,const void*))
{
int ii,jj;
double complex t = r;
double phase;
char JOBVL ='N';   // Compute Right eigenvectors
char JOBVR ='V';   // Do not compute Left eigenvectors
double complex VL[1];
int LDVL = 1; 
int LDVR = N;
int LWORK = 4*N; 
int RRWORK = 2*N;
double complex *T = (double complex*) malloc((N*N+1)*sizeof(double complex));
double complex *TT = (double complex*) malloc((N*N+1)*sizeof(double complex));
double complex *WORK =  (double complex*)malloc(LWORK*sizeof(double complex));
double complex *RWORK = (double complex*)malloc(RRWORK*sizeof(double complex));
int INFO;

for(ii=0;ii<N*N;ii++)
{
	T[ii] = S[ii]*cpow(t,z[ii]);
}

//Take the transpose of matrix T
for(ii=0;ii<N;ii++)
{
	for(jj=0;jj<N;jj++) 
	{
		TT[ii+N*jj] = T[ii*N+jj];
	}
}

//Use LAPACK to calculate all eigenvalues
zgeev_( &JOBVL, &JOBVR, &N, TT ,  &N , eigenvaluesW ,  VL, &LDVL,  eigenvectorsVR, &LDVR,   WORK,   &LWORK, RWORK, &INFO );
//Sort the eigenvalues (they're not always sorted correctly!)
qsort(eigenvaluesW,N,sizeof(eigenvaluesW[0]),fp);
 
/*
for(ii=0;ii<N;ii++)
{
	printf("%f %f: %f\n",creal(eigenvaluesW[ii]),cimag(eigenvaluesW[ii]),cabs(eigenvaluesW[ii]));
}
*/

phase = 1/M_PI*acos((creal(eigenvaluesW[0])*creal(eigenvaluesW[1])+cimag(eigenvaluesW[0])*cimag(eigenvaluesW[1]))/(cabs(eigenvaluesW[0])*cabs(eigenvaluesW[1])));
printf("Phase = %.16g\n",phase);

free(WORK);
free(RWORK);
free(T);
free(TT);

return fabs((cabs(eigenvaluesW[0])-cabs(eigenvaluesW[1]))/cabs(eigenvaluesW[1]));  
}

int compfn(const void* x, const void* y)
{
  complex double xx = *(complex double*)x;
  complex double yy = *(complex double*)y;
  if (cabs(xx) > cabs(yy)) return -1;
  if (cabs(xx) < cabs(yy)) return  1;
  return 0;
}

int main()
{
unsigned long ii, jj;
int tmp;
const int L = 6;
double r = -0.8;
int N;
unsigned long Fibon[37] = {0,1,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657,46368,75025,121393,196418,317811,514229,832040,1346269,2178309,3524578,5702887,9227465,14930352,24157817,39088169}; //valid up to L=36. L=0 is chosen as zero
int row, column, exp;
int *S;
int *z;
double complex *TransferM;
double complex *eigenVectors;
double complex *eigenValues;
char filename[35];
char widthstring[3];
FILE *fid;
time_t tic;
time_t toc;
double totaltime;
int fileinputcheck=1;
double diff;
int (*fp)(const void*,const void*) = compfn;

N = Fibon[L];

S = (int*)malloc((N*N+1) * sizeof(int));
z = (int*)malloc((N*N+1) * sizeof(int));
TransferM = (double complex*)malloc((N*N+1) * sizeof(double complex));
eigenVectors = (double complex*)malloc((N*N+1) * sizeof(double complex));
eigenValues = (double complex*)malloc((N+1) * sizeof(double complex));
if((S == NULL) & (z == NULL) & (TransferM == NULL))
{
	fprintf(stderr, "No memory left for allocating with malloc for S, z, and TransferM matrices.");
	exit(EXIT_FAILURE);
}
  
strcpy(filename,"arpack/hard_squares_fb_");
sprintf(widthstring,"%d",L);
strcat(filename,widthstring);
strcat(filename,".txt");
fid = fopen(filename,"r");
if (fid == NULL)
{
	printf("Could not open file %s\n",filename);
	exit(EXIT_FAILURE);
}
	
//Zero the matrices
for(ii=0;ii<N*N+1;ii++)
{
	S[ii]=0;
	z[ii]=0;
}

time(&tic);
//Read the transfer matrix file
while((fileinputcheck != 0) & (fileinputcheck != EOF))
{
	fileinputcheck = fscanf(fid,"%d %d %d",&row,&column,&exp);
	z[N*(row-1)+(column-1)]=exp;
	S[N*(row-1)+(column-1)]=1;
}
fclose(fid);
time(&toc);
totaltime = difftime(toc,tic);
printf("The time to read data was %gs.\n",totaltime);

time(&tic);
printf("L = %d:\n",L);
printf("The x location is %.16g\n",r);
diff = EigenRatio( eigenVectors, eigenValues,r,S,z,N,fp);
time(&toc);
totaltime = difftime(toc,tic);
printf("The time to calculate eigenvalues was %gs.\n",totaltime);
//printf("\nRelative difference = %g\n",diff);
	
return 0;
}