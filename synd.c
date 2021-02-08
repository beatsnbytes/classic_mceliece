/*
  This file is for syndrome computation
*/

#include "synd.h"

#include "params.h"
#include "root.h"
#include <sys/time.h>
#include <stdio.h>

/* input: Goppa polynomial f, support L, received word r */
/* output: out, the syndrome of length 2t */
void synd(gf *out, gf *f, gf *L, unsigned char *r)
{
	int i, j;
	gf e_inv, c;

	gf e_mat[SYS_N];

	for (j = 0; j < 2*SYS_T; j++)
		out[j] = 0;

	struct timeval start, end;
    	gettimeofday(&start, NULL);		

	eval(f, L, e_mat);

	gettimeofday(&end, 0);
    	long seconds = end.tv_sec - start.tv_sec;
    	long microseconds = end.tv_usec - start.tv_usec;
    	double elapsed_eval = seconds + microseconds*0.000001;
	sum_eval += elapsed_eval;
	times_eval = times_eval + 1;

	for (i = 0; i < SYS_N; i++)
	{
		c = (r[i/8] >> (i%8)) & 1;


		e_inv = gf_inv(gf_mul(e_mat[i],e_mat[i]));

		for (j = 0; j < 2*SYS_T; j++)
		{
			out[j] = gf_add(out[j], gf_mul(e_inv, c));
			e_inv = gf_mul(e_inv, L[i]);
		}
	}
}

