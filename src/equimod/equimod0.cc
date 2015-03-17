// This function uses the c++ version of brent from Numerical Recipes along with ARPACK++

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#define _USE_MATH_DEFINES //In order to use M_PI for pi
#include <math.h>
//#include <complex.h>
#include "arcomp.h"		//ARPACK
#include "ardscomp.h"	//ARPACK
#include "ardnsmat.h"	//ARPACK
#include "blas1c.h"		//ARPACK
#include "nr3_errno.h"		//Numerical Recipes defitions
#include "partra.h"
//#include "mpc.h" 


int compare (const void *a, const void *b)
{
	if ( abs(*(arcomplex<double>*)a) > abs(*(arcomplex<double>*)b)) return -1;
	if ( abs(*(arcomplex<double>*)a) == abs(*(arcomplex<double>*)b)) return 0;
	if ( abs(*(arcomplex<double>*)a) < abs(*(arcomplex<double>*)b)) return 1;
}

double EigenRatio(const int maxit,double *ratio, int *flag, double *phase, unsigned long long* msize, unsigned char**** M, arcomplex<double>  z, const double x, const double y, arcomplex<double> *valA, arcomplex<double> *eigenarray, const unsigned long long numeigs)
{
	unsigned long long n,m,p,q;
	
	//Defining t
	arcomplex<double> t(x,y);
	//std::cout << "\nComplex value z = " << z << std::endl;
	//std::cout << "Complex value t = " << t << std::endl;
	
	// Creating matrix A
	for(n=0ULL;n<msize[0];n++)
	{
		for (m=0ULL;m<msize[0];m++)
		{
			for (p=0ULL;p<M[n][m][0][0];p++)
			{
				valA[n*msize[0]+m] = arcomplex<double>(M[n][m][1][msize[1]*p+2],0.0)*pow(t,M[n][m][1][msize[1]*p])*pow(z,M[n][m][1][msize[1]*p+1]);
			}
		}
	}
	ARdsNonSymMatrix<arcomplex<double>, double> A(msize[0], valA);

	// Defining the problem: the two eigenvalues of A with largest magnitude
	ARluCompStdEig<double> dprob(numeigs, A);
	dprob.ChangeTol(1e-32); //default tol = 1e-16
	//printf("The default tol is %.32g\n", (double)dprob.GetTol());

	dprob.ChangeMaxit(maxit);

	//Find Eigenvalues
	dprob.FindEigenvalues();

	//Find ratio and phase and set flag
	if (dprob.EigenvaluesFound()&(dprob.ConvergedEigenvalues() == numeigs))
	{
		for(n=0ULL;n<numeigs;n++)
		{
			eigenarray[n]=dprob.Eigenvalue(n);
		}
		qsort(eigenarray,numeigs,sizeof(arcomplex<double>),compare);

		//std::cout << "E1 = " << eigenarray[0] << std::endl;
		//std::cout << "E2 = " << eigenarray[1] << std::endl;
		//std::cout << "E3 = " << eigenarray[2] << std::endl;
		*ratio = fabs(abs(eigenarray[0])-abs(eigenarray[1]))/abs(eigenarray[1]);
		*phase = acos((real(eigenarray[0])*real(eigenarray[1])+imag(eigenarray[0])*imag(eigenarray[1]))/(abs(eigenarray[0])*abs(eigenarray[1])))/M_PI;
		*flag = 0;
	}
	else 
	{
		*flag = 1;
	}

	return *ratio;
}


struct Bracketmethod {
	Doub ax,bx,cx,fa,fb,fc;
	template <class T>
	void bracket(const Doub a, const Doub b, T &func)
	{
		const Doub GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
		ax=a; bx=b;
		Doub fu;
		fa=func(ax);
		fb=func(bx);
		if (fb > fa) {
			SWAP(ax,bx);
			SWAP(fb,fa);
		}
		cx=bx+GOLD*(bx-ax);
		fc=func(cx);
		while (fb > fc) {
			Doub r=(bx-ax)*(fb-fc);
			Doub q=(bx-cx)*(fb-fa);
			Doub u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*SIGN(MAX(abs(q-r),TINY),q-r));
			Doub ulim=bx+GLIMIT*(cx-bx);
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
	inline void shft2(Doub &a, Doub &b, const Doub c)
	{
		a=b;
		b=c;
	}
	inline void shft3(Doub &a, Doub &b, Doub &c, const Doub d)
	{
		a=b;
		b=c;
		c=d;
	}
	inline void mov3(Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
		const Doub f)
	{
		a=d; b=e; c=f;
	}
};


struct Brent : Bracketmethod {
	Doub xmin,fmin;
	const Doub tol;
	Brent(const Doub toll=3.0e-8) : tol(toll) {}
	Doub minimize(double (*func)(const int,double*, int*, double*, unsigned long long*, unsigned char****, arcomplex<double>, const double, const double, arcomplex<double>*, arcomplex<double>*, const unsigned long long),const int maxit, double* ratiop, int* flagp, double* phasep, unsigned long long* msize, unsigned char**** M, arcomplex<double> z, const double angle, arcomplex<double> *valA, arcomplex<double>* eigenarray, const unsigned long long numeigs)
	{
		const Int ITMAX=100;
		const Doub CGOLD=0.3819660;
		const Doub ZEPS=numeric_limits<Doub>::epsilon()*1.0e-3;
		Doub a,b,d=0.0,etemp,fu,fv,fw,fx;
		Doub p,q,r,tol1,tol2,u,v,w,x,xm;
		Doub e=0.0;
		Doub fval;
	
		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,z,angle,x,valA,eigenarray,numeigs);
		if (*flagp==0)
		{
			fw=fv=fx=fval;
		}
		else
		{
			//std::cout << "ARPACK did not converge" << std::endl;
			return fval;
		}
		for (Int iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*abs(x)+ZEPS);
			if (abs(x-xm) <= (tol2-0.5*(b-a))) {
				fmin=fx;
				return xmin=x;
			}
			if (abs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=abs(q);
				etemp=e;
				e=d;
				if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x)
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
			u=(abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fval = (*func)(maxit,ratiop,flagp,phasep,msize,M,z,angle,u,valA,eigenarray,numeigs);
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
		cout << "Too many iterations in brent" << std::endl;
		return fval;
	}
};




int main()
{
	unsigned char N=4;
	unsigned long long int Rtotal = 100000; 	//Total number of points
	double ratiotol=1e-8;	//Tolerance for Convergence in Brent
	int maxit=50000;	//Maximum number of ARPACK iterations (default used was 50000)
	double xmin = 0;
	double xmax = 4;
	double ymin = 0;
	double ymax = 4;
	unsigned long long numeigs=11;
	
	//power law transformation of uniform distribution of random numbers. 
	//flatness of 1 is flat. flatness of n means range*(rand^n), where rand is in the interval [0,1].
	double xflatness = 2;	
	double yflatness = 2;
	
	unsigned long long ii, jj, tmp;	// Counters and tmp
	unsigned long long arerrors=0, brenterrors=0, localmin=0;

	unsigned long long msize[2];      	// TM size
	partra_matrix M;
	arcomplex<double> z(0.99,0.0);
	double x,y1,y2,y3;
	double r;
	arcomplex<double>* valA;   	// pointer to an array that stores values of A
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

	double (*objtfnpt)(const int,double*, int*, double*, unsigned long long*, unsigned char****, arcomplex<double>, const double, const double, arcomplex<double>*,arcomplex<double>*,const unsigned long long) = &EigenRatio;

	std::cout << "Choose an output file number: " ;
	std::cin >> n;


	srand (time(NULL)); //Seed the random number generator

	//Create transfer matrix
	flag = if_sq_c_f(&M,msize,filename1,N);
	if (flag!=0)
	{
		return 0;
	}
	
	//Create array for eigenvalues
	if(msize[0]<11)
	{
		numeigs = msize[0]-1;
	}
	arcomplex<double> *eigenarray = new arcomplex<double>[numeigs];

	//Open the output file
	sprintf(filename2,"equimod_%s_%.4f_%llu.txt",filename1,abs(z),n);
	fid = fopen(filename2,"w");
	if (fid == NULL)
	{
		printf("ERROR: Could not open output file %s\n",filename2);
		matrix_free(M,msize);
		return 0;
	}
	setvbuf(fid,NULL,_IOLBF,32); 


	// Defining Matrix valA
	valA  = new arcomplex<double>[msize[0]*msize[0]];
	
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
			r = brent.minimize(objtfnpt,maxit,&ratio,&flag,&phase,msize,M,z,x,valA,eigenarray,numeigs);

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
				arerrors++;
			}
			else if (flag == 2)
			{
				brenterrors++;
			}
	}
	time(&toc);
	totaltime = difftime(toc,tic);

	printf("\n\nThe width N = %d\n",N);
	printf("ARPACK errors: %llu/%llu\n",arerrors,Rtotal);
	printf("brent errors: %llu/%llu\n",brenterrors,Rtotal);
	printf("local min errors: %llu/%llu\n",localmin,Rtotal);
	printf("TOTAL CONVERGED: %llu/%llu\n",Rtotal-arerrors-brenterrors-localmin,Rtotal);
	std::cout << "The total time to calculate equimod was " << totaltime << "s" << std::endl;
	std::cout << std::endl;

	fclose(fid);
	delete[] valA;
	matrix_free(M,msize);
	delete[] eigenarray;
	return 0;
}

