// This function uses the c++ version of brent from Numerical Recipes along with ARPACK++

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#define _USE_MATH_DEFINES //In order to use M_PI for pi
#include <math.h>
#include <complex.h>
//#include <limit.hs>
#include "partra.h"

//Numerical Recipes:
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int compare (const void *a, const void *b)
{
	if ( cabs(*(double complex*)a) > cabs(*(double complex*)b)) return -1;
	if ( cabs(*(double complex*)a) == cabs(*(double complex*)b)) return 0;
	if ( cabs(*(double complex*)a) < cabs(*(double complex*)b)) return 1;
}

double EigenRatio(const int maxit,double *ratio, int *flag, double *phase, unsigned long long* msize, double**** M, const double row, const double x, double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *T, double complex *TT, double complex *WORK, double complex *RWORK)
{
  unsigned long long n,m,p,q;
  double complex t = row+I*x;
  char JOBVL ='N';   // Compute Right eigenvectors
  char JOBVR ='V';   // Do not compute Left eigenvectors
  double complex VL[1];
  int s = msize[0];
  int LDVL = 1; 
  int LDVR = msize[0];
  int LWORK = 4*msize[0]; 
  int INFO;
  
 
  for(n=0ULL;n<msize[0];n++)
  {
    for (m=0ULL;m<msize[0];m++)
    {//printf("%f %f\n",M[n][m][1][0],M[n][m][1][1]);
      for (p=0ULL;p<M[n][m][0][0];p++)
      {
	T[n*msize[0]+m] =T[n*msize[0]+m] + M[n][m][1][msize[1]*p+1]*cpow(t,M[n][m][1][msize[1]*p]);
	//printf("T[%llu]=%f+I*%f\n",n*msize[0]+m,crealf(T[n*msize[0]+m]),cimagf(T[n*msize[0]+m]));
      }
    }
  }
  //printf("t = %f+I*%f, M[0][2] = %f, T[0] = %f+I*%f\n",crealf(t),cimagf(t),M[0][2][1][1],crealf(cpowf(t,2)),cimagf(cpowf(t,2)));
  
  //Take the transpose of matrix T
  for(n=0;n<msize[0];n++)
  {
    for(m=0;m<msize[0];m++) 
    {
      TT[n+msize[0]*m] = T[n*msize[0]+m];
    }
  }

  //Use LAPACK to calculate all eigenvalues
  zgeev_( &JOBVL, &JOBVR, &s, TT ,  &s , eigenvaluesW ,  VL, &LDVL,  eigenvectorsVR, &LDVR,   WORK,   &LWORK, RWORK, &INFO );
  if (INFO==0)
  {
    qsort(eigenvaluesW,msize[0],sizeof(double complex),compare);
    *phase = 1/M_PI*acos((creal(eigenvaluesW[0])*creal(eigenvaluesW[1])+cimag(eigenvaluesW[0])*cimag(eigenvaluesW[1]))/(cabs(eigenvaluesW[0])*cabs(eigenvaluesW[1])));
    *ratio = fabs((cabs(eigenvaluesW[0])-cabs(eigenvaluesW[1]))/cabs(eigenvaluesW[1]));
    *flag = 0;
    printf("\nx = %f, y = %f, e1=%f, e2=%f, ratio = %f",row,x,cabs(eigenvaluesW[0]),cabs(eigenvaluesW[1]),*ratio);
  }
  else 
  {
    *flag = 1;
  }  

  return *ratio;  
}

double brent(double ax,double bx,double cx,double tol,double (*func)(const int,double*, int*, double*, unsigned long long*, double****, const double, const double, double complex *,double complex *,double complex *,double complex *, double complex *, double complex *),const int maxit, double* ratiop, int* flagp, double* phasep, unsigned long long* msize, double**** M, const double row, double complex *eigenvectors, double complex *eigenvalues, double complex *T,double complex *TT,double complex *WORK,double complex *RWORK)
{
      int iter;
      double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
      double e=0.0;
      double xmin,fmin,fval;

      a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,row,x,eigenvectors,eigenvalues,T,TT,WORK,RWORK);
	if (*flagp==0)
	{
	  fw=fv=fx=fval;
	}
	else
	{
	  //std::cout << "ARPACK did not converge" << std::endl;
	  return fval;
	}
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			return xmin=x;;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,row,u,eigenvectors,eigenvalues,T,TT,WORK,RWORK);
		if (*flagp==0)
		{
		    fu=fval;
		}
		else
		{
		    //std::cout << "ARPACK did not converge" << std::endl;
		    return fval;
		}
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	*flagp = 2;
	printf("Too many iterations in brent\n");
	return fval;
}




int main()
{//TC0+:4:4,5:6,6:11,7:16,8:28,9:44,10:76,11:124,12:222,13:378
  //TF+:4:9,5:18,6:34,7:70,8:134,9:270,10:526
  // TF:4:14,5:30,6:62,7:126,8:254,9:510,10:1022
unsigned char N=4;
unsigned long long int Rtotal = 1; 	//Total number of points
double ratiotol=1e-14;	//Tolerance for Convergence in Brent. Default = 3e-8
int maxit=50000;	//Maximum number of ARPACK iterations (default used was 50000)
double xmin = 0.;
double xmax = 3.;
double ymin = 0;
double ymax = 2.0;
double z = 1.;		//Field H value. H=0 is z=1. H=infinity is z=0

//power law transformation of uniform distribution of random numbers. 
//flatness of 1 is flat. flatness of n means range*(rand^n), where rand is in the interval [0,1].
double xflatness = 1.05;	
double yflatness = 1.05;

unsigned long long ii, jj, tmp;	// Counters and tmp
unsigned long long laperrors=0, brenterrors=0, localmin=0;

unsigned long long misize[2];      	// TM size
unsigned long long msize[2];
partra_matrix Mi;
partra_matrix_d M;
double x,y1,y2,y3;
double r;
double ratio;
double phase;
int flag;

unsigned long long n;
int scheck;
char filename1[256];
char filename2[256];
char widthstring[3];
FILE *fid;
int fileinputcheck=1;

time_t tic;
time_t toc;
double totaltime;


double (*objtfnpt)(const int,double*, int*, double*, unsigned long long*, double****, const double, const double,double complex *,double complex *,double complex *,double complex *,double complex *,double complex *) = &EigenRatio;

//Output file name
printf("Choose an output file number: ");
scheck=scanf("%llu",&n);
if (scheck<=0)
{
  printf("ERROR: %s",strerror(errno));
  return 0;
}

srand (time(NULL)); //Seed the random number generator

//Create transfer matrix
flag = if_sq_f_r(&Mi,misize,filename1,N);
if (flag!=0)
{
  return 0;
}

//Substitute value
char subs[256]="x";
flag = matrix_sub_d(&M,msize,Mi,misize,subs,z);
if (flag!=0)
{
  matrix_free(Mi,misize);
  return 0;
}
matrix_free(Mi,misize);

//Create LAPACK vectors
double complex *eigenvectors = (double complex*)calloc((msize[0]*msize[0]+1),sizeof(double complex));
double complex *eigenvalues = (double complex*)calloc((msize[0]+1),sizeof(double complex));
double complex *T = (double complex*) calloc((msize[0]*msize[0]+1),sizeof(double complex));
double complex *TT = (double complex*) calloc((msize[0]*msize[0]+1),sizeof(double complex));
double complex *WORK =  (double complex*)calloc(4*msize[0],sizeof(double complex));
double complex *RWORK = (double complex*)calloc(2*msize[0],sizeof(double complex));

//Open the output file
sprintf(filename2,"equimod_%s_%.4f_%llu.txt",filename1,fabs(z),n);
fid = fopen(filename2,"w");
if (fid == NULL)
{
  printf("ERROR: Could not open output file %s\n",filename2);
  matrix_free_d(M,msize);
  return 0;
}
setvbuf(fid,NULL,_IOLBF,32); 



//Main Loop
time(&tic);
printf("The width is %hhu\n",N);
for (ii=0;ii<Rtotal;ii++)
{
  printf("\n%llu",ii);
  x = (xmax-xmin)*pow(rand()/((double)RAND_MAX),xflatness)+xmin;
  y1 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;
  y2 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;
  y3 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;

  // Finding eigenvalues
  r = brent(y1,y2,y3,ratiotol,objtfnpt,maxit,&ratio,&flag,&phase,msize,M,x,eigenvectors,eigenvalues,T,TT,WORK,RWORK);

  // Printing solution.
  //printf(" flag = %d, r = %f, ratio = %f",flag,r, ratio);
  if ((flag == 0) &(ratio<ratiotol))
  {
    fprintf(fid,"%.8g %.8g %.8g\n",x,r,phase);
    printf("*");
  }
  else if ((flag == 0) &(ratio>ratiotol))
  {
    localmin++;
  }
  else if (flag == 1)
  {
    laperrors++;
  }
  else if (flag == 2)
  {
    brenterrors++;
  }
}
time(&toc);
totaltime = difftime(toc,tic);

printf("\n\nThe width N = %d\n",N);
printf("LAPACK errors: %llu/%llu\n",laperrors,Rtotal);
printf("brent errors: %llu/%llu\n",brenterrors,Rtotal);
printf("local min errors: %llu/%llu\n",localmin,Rtotal);
printf("TOTAL CONVERGED: %llu/%llu\n",Rtotal-laperrors-brenterrors-localmin,Rtotal);
printf("The total time to calculate equimod was %fs\n\n",totaltime);

fclose(fid);
matrix_free_d(M,msize);
free(eigenvectors);
free(eigenvalues);
free(WORK);
free(RWORK);
free(T);
free(TT);
return 0;
}

