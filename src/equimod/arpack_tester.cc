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
//#include "nr3.h"		//Numerical Recipes defitions
//#include "mpc.h" 
#include "partra.h"

int compare (const void *a, const void *b)
{
	if ( abs(*(arcomplex<double>*)a) > abs(*(arcomplex<double>*)b)) return -1;
	if ( abs(*(arcomplex<double>*)a) == abs(*(arcomplex<double>*)b)) return 0;
	if ( abs(*(arcomplex<double>*)a) < abs(*(arcomplex<double>*)b)) return 1;
}


double EigenRatio(const int maxit,double *ratio, int *flag, double *phase, unsigned long long* msize, unsigned char**** M, arcomplex<double>  z, const double x, const double y, arcomplex<double> *valA)
{
	unsigned long long n,m,p;
	unsigned char numeigs=11;
	if(msize[0]<11)
	{
		numeigs = msize[0]-1;
	}

	arcomplex<double> *eigenarray = new arcomplex<double>[numeigs];
	
	//Defining t
	arcomplex<double> t(x,y);
	std::cout << "\nComplex value z = " << z << std::endl;
	std::cout << "Complex value t = " << t << std::endl;
	
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

	delete[] eigenarray;
	return *ratio;
}



int main()
{
	unsigned char N=4;
	unsigned long long int Rtotal = 1; 	//Total number of points
	double ratiotol=1e-8;	//Tolerance for Convergence in Brent
	int maxit=50000;	//Maximum number of ARPACK iterations
	double xmin = -1.01;
	double xmax = -0.99;
	double ymin = 0;
	double ymax = 0.01;
	
	//power law transformation of uniform distribution of random numbers. 
	//flatness of 1 is flat. flatness of n is rand^n.
	double xflatness = 1;	
	double yflatness = 1;
	
	unsigned long long ii, jj, tmp;	// Counters and tmp
	unsigned long long arerrors=0, brenterrors=0, localmin=0;

	unsigned long long msize[2];      	// TM size
	partra_matrix M;
	arcomplex<double> z(0.2,0.4);


	double x,y;
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

	//std::cout << "Choose a width L: " ;
	//std::cin >> L;
	std::cout << "Choose an output file number: " ;
	std::cin >> n;


	srand (time(NULL)); //Seed the random number generator

	//Create transfer matrix
	flag = if_sq_c_f(&M,msize,filename1,N);
	if (flag!=0)
	{
		return 0;
	}

	//Open the output file
	sprintf(filename2,"equimod_%s_%llu.txt",filename1,n);
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
	
	
	//Main Loop
	time(&tic);
	printf("The width is %hhu\n",N);
	for (ii=0;ii<Rtotal;ii++)
	{std::cout<<std::endl<<ii;
			x = (xmax-xmin)*pow(rand()/((double)RAND_MAX),xflatness)+xmin;
			y = (ymax-ymin)*pow(rand()/((double)RAND_MAX),yflatness)+ymin;
			
			// Finding eigenvalues
			EigenRatio(maxit,&ratio,&flag,&phase,msize,M,z,x,y,valA);
			
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
	return 0;
}

