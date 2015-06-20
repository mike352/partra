#include "partra.h"
#include <complex.h>

int main(void)
{
unsigned char N=4;
char dirname[256] = "data";
unsigned char flag;

partra_matrix M;
unsigned long long msize[2],mosize[2];
char filename1[256],filename2[256];

//Create and output first file
flag = if_sq_f_r_s_file(N,dirname); 
if (flag!=0)
{
	return 0;
}

//Second file - outputs should be equal
flag = if_sq_f_r_s(&M,msize,filename1,N);
sprintf(filename2,"%s/_%s.txt",dirname,filename1);
flag = matrix_fprintf(M,msize,filename2);

//Substitution tests
partra_matrix_d Mo1;
char subs1[256]="sdx";
double z=0.23;
flag = matrix_sub_d(&Mo1,mosize,M,msize,subs1,z);
if (flag!=0)
{
  matrix_free(M,msize);
  return 0;
}
sprintf(filename2,"%s/%s_x23.txt",dirname,filename1);
flag = matrix_fprintf_d(Mo1,mosize,filename2,16);
matrix_free_d(Mo1,mosize);


partra_matrix_d Mo2;
char subs2[256]="sdu";
z=0.32;
flag = matrix_sub_d(&Mo2,mosize,M,msize,subs2,z);
if (flag!=0)
{
  matrix_free(M,msize);
  return 0;
}
sprintf(filename2,"%s/%s_u32.txt",dirname,filename1);
flag = matrix_fprintf_d(Mo2,mosize,filename2,16);

partra_matrix_d Mo22;
unsigned long long mo2size[2];
z=0.32;
char subs22[256]="s";
flag = matrix_sub_d_d(&Mo22,mo2size,Mo2,mosize,subs22,z);
if (flag!=0)
{
  matrix_free(M,msize);
  matrix_free_d(Mo2,mosize);
  return 0;
}
sprintf(filename2,"%s/%s_u32x32.txt",dirname,filename1);
flag = matrix_fprintf_d(Mo22,mo2size,filename2,16);
matrix_free_d(Mo2,mosize);
matrix_free_d(Mo22,mosize);


partra_matrix_d Mo3;
char subs3[256]="sdux";
double z1=0.23,z2=0.32;
flag = matrix_sub_d(&Mo3,mosize,M,msize,subs3,z1,z2);
if (flag!=0)
{
  matrix_free(M,msize);
  return 0;
}
sprintf(filename2,"%s/%s_u23x32.txt",dirname,filename1);
flag = matrix_fprintf_d(Mo3,mosize,filename2,16);
matrix_free_d(Mo3,mosize);

partra_matrix_d Mo4;
char subs4[256]="sdxu";
z1=0.23,z2=0.32;
flag = matrix_sub_d(&Mo4,mosize,M,msize,subs4,z1,z2);
if (flag!=0)
{
  matrix_free(M,msize);
  return 0;
}
sprintf(filename2,"%s/%s_u32x23.txt",dirname,filename1);
flag = matrix_fprintf_d(Mo4,mosize,filename2,16);
matrix_free_d(Mo4,mosize);

/*
fid = fopen(filename2,"w");
for (n=0ULL;n<msize[0];n++)
{
	for (m=0ULL;m<msize[0];m++)
	{
		fprintf(fid,"%hhu\n",M[n][m][1][0]);
	}
}
fclose(fid);
*/


matrix_free(M,msize);
return 0;
}