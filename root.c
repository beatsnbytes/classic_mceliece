/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include "root.h"
#include "params.h"
#include "gf.h"

#include <stdio.h>
#include <sys/time.h>

double sum_eval=0;
int times_eval=0;

/* input: polynomial f and field element a */
/* return f(a) */
void eval(gf *f, gf *a, gf* out)
{

	int i, j, k;
	gf r, r_tmp;
	uint32_t tmp;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;
	
	for (k = 0; k < SYS_N; k++)
	{

		r = f[ SYS_T ];
		for (i = SYS_T-1; i >= 0; i--)
		{

			t0 = r;
			t1 = a[k];

			tmp = t0 * (t1 & 1);

			for (j = 1; j < GFBITS; j++)
				tmp ^= (t0 * (t1 & (1 << j)));

			t = tmp & 0x7FC000;;
			tmp ^= (t >> 9) ^ (t >> 12);

			t = tmp & 0x3000;
			tmp ^= (t >> 9) ^ (t >> 12);

			r_tmp = tmp & ((1 << GFBITS)-1);

			r = r_tmp ^ f[i];
		}

		out[k]=r;
	}

	//return r;
}

/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root(gf *out, gf *f, gf *L)
{


	struct timeval start, end;
    gettimeofday(&start, NULL);

	eval(f, L, out);

	gettimeofday(&end, 0);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed_eval = seconds + microseconds*0.000001;
	sum_eval += elapsed_eval;
	times_eval = times_eval + 1;
}

