/*
  This file is for Niederreiter encryption
*/

#include "encrypt.h"

#include "util.h"
#include "params.h"
#include "randombytes.h"

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#include "gf.h"

double sum_syndrome = 0;
int times_syndrome = 0;

static inline unsigned char same_mask(uint16_t x, uint16_t y)
{
	uint32_t mask;

	mask = x ^ y;
	mask -= 1;
	mask >>= 31;
	mask = -mask;

	return mask & 0xFF;
}

/* output: e, an error vector of weight t */
static void gen_e(unsigned char *e)
{
	int i, j, eq;

	uint16_t ind[ SYS_T ];
	unsigned char bytes[ sizeof(ind) ];
	unsigned char mask;	
	unsigned char val[ SYS_T ];	

	while (1)
	{
		randombytes(bytes, sizeof(bytes));

		for (i = 0; i < SYS_T; i++)
			ind[i] = load_gf(bytes + i*2);

		// check for repetition

		eq = 0;

		for (i = 1; i < SYS_T; i++) 
			for (j = 0; j < i; j++)
				if (ind[i] == ind[j]) 
					eq = 1;

		if (eq == 0)
			break;
	}

	for (j = 0; j < SYS_T; j++)
		val[j] = 1 << (ind[j] & 7);

	for (i = 0; i < SYS_N/8; i++) 
	{
		e[i] = 0;

		for (j = 0; j < SYS_T; j++)
		{
			mask = same_mask(i, (ind[j] >> 3));

			e[i] |= val[j] & mask;
		}
	}
}

/* input: public key pk, error vector e */
/* output: syndrome s */
static void syndrome(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	unsigned char b, row[SYS_N/8];
	const unsigned char *pk_ptr = pk;

	int i, j;

	for (i = 0; i < SYND_BYTES; i++)
		s[i] = 0;

	for (i = 0; i < PK_NROWS; i++)	
	{
		for (j = 0; j < SYS_N/8; j++) 
			row[j] = 0;

		for (j = 0; j < PK_ROW_BYTES; j++) 
			row[ SYS_N/8 - PK_ROW_BYTES + j ] = pk_ptr[j];

		row[i/8] |= 1 << (i%8);
		
		b = 0;
		for (j = 0; j < SYS_N/8; j++)
			b ^= row[j] & e[j];

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

		s[ i/8 ] |= (b << (i%8));

		pk_ptr += PK_ROW_BYTES;
	}
}

void encrypt(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	gen_e(e);

#ifdef KAT
  {
    int k;
    printf("encrypt e: positions");
    for (k = 0;k < SYS_N;++k)
      if (e[k/8] & (1 << (k&7)))
        printf(" %d",k);
    printf("\n");
  }
#endif

	struct timeval start_syndrome, end_syndrome;
    gettimeofday(&start_syndrome, NULL);

	syndrome(s, pk, e);

	gettimeofday(&end_syndrome, 0);
    long seconds = end_syndrome.tv_sec - start_syndrome.tv_sec;
    long microseconds = end_syndrome.tv_usec - start_syndrome.tv_usec;
    double elapsed_syndrome = seconds + microseconds*0.000001;
	sum_syndrome += elapsed_syndrome;
	times_syndrome = times_syndrome + 1;
}

