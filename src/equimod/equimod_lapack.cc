// This function uses the c++ version of brent from Numerical Recipes along with ARPACK++

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#define _USE_MATH_DEFINES //In order to use M_PI for pi
#include <math.h>
#include <complex.h>
#include <limits>
#include "partra.h"
#include <iostream> //C++
#include "gmp.h"
//#include "mpc.h" 


int compare (const void *a, const void *b)
{
	if ( cabs(*(double complex*)a) > cabs(*(double complex*)b)) return -1;
	if ( cabs(*(double complex*)a) == cabs(*(double complex*)b)) return 0;
	if ( cabs(*(double complex*)a) < cabs(*(double complex*)b)) return 1;
}

double EigenRatio(const int maxit,double *ratio, int *flag, double *phase, unsigned long long* msize, double**** M, const double x, const double y, double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *T, double complex *TT, double complex *WORK, double complex *RWORK, const unsigned long long numeigs)
{
  unsigned long long n,m,p,q;
  double complex t = x+I*y;
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
	T[n*msize[0]+m] =T[n*msize[0]+m] + (M[n][m][1][msize[1]*p+1],0.0)*cpow(t,M[n][m][1][msize[1]*p]);
      }
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
  }
  else 
  {
    *flag = 1;
  }  

  return *ratio;  
}



static const double NaN = std::numeric_limits<double>::quiet_NaN();

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}


struct Bracketmethod {
  double ax,bx,cx,fa,fb,fc;
  template <class T>
  void bracket(const double a, const double b, T &func)
  {
    const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
    ax=a; bx=b;
    double fu;
    fa=func(ax);
    fb=func(bx);
    if (fb > fa) {
      SWAP(ax,bx);
      SWAP(fb,fa);
    }
    cx=bx+GOLD*(bx-ax);
    fc=func(cx);
    while (fb > fc) {
      double r=(bx-ax)*(fb-fc);
      double q=(bx-cx)*(fb-fa);
      double u=bx-((bx-cx)*q-(bx-ax)*r)/
	      (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
      double ulim=bx+GLIMIT*(cx-bx);
      if ((bx-u)*(u-cx) > 0.0) {
	fu=func(u);
	if (fu < fc) {
	  ax=bx;
	  bx=u;
	  fa=fb;
	  fb=fu;
	  return;
	} else if (fu > fb) {
	  cx=u;
	  fc=fu;
	  return;
	}
	u=cx+GOLD*(cx-bx);
	fu=func(u);
      } else if ((cx-u)*(u-ulim) > 0.0) {
	fu=func(u);
	if (fu < fc) {
	  shft3(bx,cx,u,u+GOLD*(u-cx));
	  shft3(fb,fc,fu,func(u));
	}
      } else if ((u-ulim)*(ulim-cx) >= 0.0) {
	u=ulim;
	fu=func(u);
      } else {
	u=cx+GOLD*(cx-bx);
	fu=func(u);
      }
      shft3(ax,bx,cx,u);
      shft3(fa,fb,fc,fu);
    }
  }
  inline void shft2(double &a, double &b, const double c)
  {
    a=b;
    b=c;
  }
  inline void shft3(double &a, double &b, double &c, const double d)
  {
    a=b;
    b=c;
    c=d;
  }
  inline void mov3(double &a, double &b, double &c, const double d, const double e,
	  const double f)
  {
    a=d; b=e; c=f;
  }
};


struct Brent : Bracketmethod {
  double xmin,fmin;
  const double tol;
  Brent(const double toll=3.0e-8) : tol(toll) {}
  double minimize(double (*func)(const int,double*, int*, double*, unsigned long long*, double****, const double, const double, double complex *,double complex *,double complex *,double complex *, double complex *, double complex *, const unsigned long long),const int maxit, double* ratiop, int* flagp, double* phasep, unsigned long long* msize, double**** M, const double angle, double complex *eigenvectors, double complex *eigenvalues, double complex *T,double complex *TT,double complex *WORK,double complex *RWORK, const unsigned long long numeigs)
  {
      const int ITMAX=100;
      const double CGOLD=0.3819660;
      const double ZEPS=std::numeric_limits<double>::epsilon()*1.0e-3;
      double a,b,d=0.0,etemp,fu,fv,fw,fx;
      double p,q,r,tol1,tol2,u,v,w,x,xm;
      double e=0.0;
      double fval;

      a=(ax < cx ? ax : cx);
      b=(ax > cx ? ax : cx);
      x=w=v=bx;
      fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,angle,x,eigenvectors,eigenvalues,T,TT,WORK,RWORK,numeigs);
      if (*flagp==0)
      {
	fw=fv=fx=fval;
      }
      else
      {
	//std::cout << "ARPACK did not converge" << std::endl;
	return fval;
      }
      for (int iter=0;iter<ITMAX;iter++) {
	xm=0.5*(a+b);
	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	  fmin=fx;
	  return xmin=x;
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
	  if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x)
			  || p >= q*(b-x))
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
	fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,angle,u,eigenvectors,eigenvalues,T,TT,WORK,RWORK,numeigs);
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
	    shft3(v,w,x,u);
	    shft3(fv,fw,fx,fu);
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
      std::cout << "Too many iterations in brent" << std::endl;
      return fval;
  }
};




int main()
{//TX0+:4:4,5:6,6:11,7:16,8:28,9:44,10:76,11:124,12:222,13:378
  //TF+:4:9,5:18,6:34,7:70,8:134,9:270,10:526
  // TF:4:14,5:30,6:62,7:126,8:254,9:510,10:1022
unsigned char N=4;
unsigned long long int Rtotal = 500000; 	//Total number of points
double ratiotol=1e-14;	//Tolerance for Convergence in Brent
int maxit=50000;	//Maximum number of ARPACK iterations (default used was 50000)
double xmin = -3.;
double xmax = 3.;
double ymin = 0;
double ymax = 4.0;
double z = 1.;		//Field H value. H=0 is z=1. H=infinity is z=0
unsigned long long numeigs=14;

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
char filename1[256];
char filename2[256];
char widthstring[3];
FILE *fid;
int fileinputcheck=1;

time_t tic;
time_t toc;
double totaltime;


double (*objtfnpt)(const int,double*, int*, double*, unsigned long long*, double****, const double, const double,double complex *,double complex *,double complex *,double complex *,double complex *,double complex *,const unsigned long long) = &EigenRatio;


std::cout << "Choose an output file number: " ;
std::cin >> n;


srand (time(NULL)); //Seed the random number generator

//Create transfer matrix
flag = if_sq_f_f(&Mi,misize,filename1,N);
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
double complex *eigenvectors = (double complex*)malloc((msize[0]*msize[0]+1) * sizeof(double complex));
double complex *eigenvalues = (double complex*)malloc((msize[0]+1) * sizeof(double complex));
double complex *T = (double complex*) malloc((msize[0]*msize[0]+1)*sizeof(double complex));
double complex *TT = (double complex*) malloc((msize[0]*msize[0]+1)*sizeof(double complex));
double complex *WORK =  (double complex*)malloc(4*msize[0]*sizeof(double complex));
double complex *RWORK = (double complex*)malloc(2*msize[0]*sizeof(double complex));

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


//Initializing optimizer brent
Brent brenttest;
Brent brent(ratiotol); //Default = 3e-8


//Main Loop
time(&tic);
printf("The width is %hhu\n",N);
for (ii=0;ii<Rtotal;ii++)
{std::cout<<std::endl<<ii;
  x = (xmax-xmin)*pow(rand()/((double)RAND_MAX),xflatness)+xmin;
  y1 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;
  y2 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;
  y3 = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;

  //Brent brent;
  brent.ax = y1;
  brent.bx = y2;
  brent.cx = y3;

  // Finding eigenvalues
  r = brent.minimize(objtfnpt,maxit,&ratio,&flag,&phase,msize,M,x,eigenvectors,eigenvalues,T,TT,WORK,RWORK,numeigs);

  // Printing solution.

  if ((flag == 0) &(ratio<ratiotol))
  {	
	  {
		  fprintf(fid,"%.8g %.8g %.8g\n",x,r,phase);
		  //printf("x=%.8g, y=%.8g\n",x,r);
		  std::cout<<"*";
	  }
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
std::cout << "The total time to calculate equimod was " << totaltime << "s" << std::endl;
std::cout << std::endl;

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

