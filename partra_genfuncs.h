#ifndef PARTRA_GENFUNCS_H
#define PARTRA_GENFUNCS_H

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>


/*This function performs a circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_single_lshift(unsigned long long number, const unsigned char width)
{
	//(shift left of number without Nth digit) + (Nth digit shifted to the front)
	number = ((~(~0ULL<<(width-1))& number)<<1) + (((1ULL<<(width-1))& number)>>(width-1));
	return number;
}

/*This function performs a binned circular shift to the left of all bits in number of a particular width, erasing any bits to the left of the given width.*/
unsigned long long circ_bin_lshift(unsigned long long number, const unsigned char N, const unsigned char bin)
{
	unsigned char n;
	for (n=0;n<bin;n++)
	{
		//(shift left of number without width digit) + (width digit shifted to the front)
		number = ((~(~0ULL<<(N*bin-1))& number)<<1) + (((1ULL<<(N*bin-1))& number)>>(N*bin-1));
	}
	return number;
}

/*This function reflects all the bits for a number of a particular width*/
unsigned long long bit_reflection(unsigned long long number, const unsigned char N)
{
	unsigned long long reflection=number;
	unsigned char s = N-1;
	//Bit twiddling from https://graphics.stanford.edu/~seander/bithacks.html
	for (number >>= 1; number; number >>= 1)
	{   
		reflection <<= 1;
		reflection |= number & 1ULL;
		s--;
	}
	reflection <<= s; 
	return (reflection & ~((~0ULL<<(N-1))<<1)); //added to algorithm to clear any set bits beyond N, done in two steps to handle case when N equals the size of long long
}

/*This function reflects all the bits for a number of a particular width and bin size bin*/
unsigned long long bit_reflection_bin(unsigned long long number, const unsigned char N, const unsigned char bin)
{
	unsigned char ii;
	unsigned long long reflection=0ULL;
	for (ii=0;ii<N;ii++)
	{
		//reflection=reflection + (select out ii-th bin and shift to starting position, shift appropriately)
		reflection = reflection + (((number & (((1ULL<<bin)-1ULL)<<bin*ii))>>bin*ii) << (bin*(N-1-ii)));
	}
	return reflection;
}

/*This function adds up all of the set bits in a number x*/
unsigned char bit_sum(unsigned long long x)
{
	unsigned char sum;
	//Standard bit twiddling
	for (sum=0;x;sum++) x &= x-1ULL;
	return sum;
}







#endif