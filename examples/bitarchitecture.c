#include <stdio.h>
#include <limits.h>

//In C, sizeof is defined as how many multiples of char a type is. 
//Therefore, we first independently determine the size of char.
//Next we independently determine the size of long longs. 
//An independent test of long doubles has not been implemented yet. 
int main(void)
{
unsigned char b=1,c=1,ll=1;
unsigned long long d=1;
for (c=1;c<256;c++)
{
	b=(b<<1);
	if (b==0) break;
}
if (c!=CHAR_BIT)
{
	printf("The C Macro CHAR_BIT is incorrect\n");
}
printf("Size of char = %d bits\n",c);


for (ll=1;ll<256;ll++)
{
	d=(d<<1);
	if (d==0) break;
}
if (ll!=c*sizeof(unsigned long long))
{
	printf("The sizeof operator does not give the correct size of long longs.\n");
}
printf("Size of long long = %d bits\n",ll);


printf("Size of long double = %lu bits\n",c*sizeof(long double));

return 0;
}