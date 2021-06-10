/*
  This file is for Niederreiter decryption
*/

#include <stdio.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include <sys/time.h>
#include "decrypt.h"

#include "params.h"
#include "encrypt.h"
#include "benes.h"
#include "util.h"
#include "synd.h"
#include "root.h"
#include "gf.h"
#include "bm.h"
#include "custom_util.h"

double sum_total_synd=0.0;
int times_total_synd=0;

/* Niederreiter decryption with the Berlekamp decoder */
/* intput: sk, secret key */
/*         c, ciphertext */
/* output: e, error vector */
/* return: 0 for success; 1 for failure */

//double sum_synd= 0.0;
//int times_synd= 0;

int decrypt(unsigned char *e, const unsigned char *sk, const unsigned char *c)
{
	int i, w = 0; 
	uint16_t check;	

	unsigned char r[ SYS_N/8 ];

	gf g[ SYS_T+1 ];
	gf L[ SYS_N ];

	gf s[ SYS_T*2 ];
	gf s_cmp[ SYS_T*2 ];
	gf locator[ SYS_T+1 ];
	gf images[ SYS_N ];

	gf t;

	//

	for (i = 0; i < SYND_BYTES; i++)       r[i] = c[i];
	for (i = SYND_BYTES; i < SYS_N/8; i++) r[i] = 0;

	for (i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; } g[ SYS_T ] = 1;

	support_gen(L, sk);

#ifdef TIME_MEASUREMENT
  	struct timeval start_synd, end_synd;
  	gettimeofday(&start_synd, NULL);
#endif
#ifdef SYND_KERNEL
	synd_host(s, g, L, r);
#endif
#ifndef SYND_KERNEL
	synd_sw_host(s, g, L, r);
#endif
#ifdef TIME_MEASUREMENT
	gettimeofday(&end_synd, NULL);
	get_event_time(&start_synd, &end_synd, &sum_total_synd, &times_total_synd);
#endif

	bm(locator, s);

	root(images, locator, L);

	//
	
	for (i = 0; i < SYS_N/8; i++) 
		e[i] = 0;

	for (i = 0; i < SYS_N; i++)
	{
		t = gf_iszero(images[i]) & 1;

		e[ i/8 ] |= t << (i%8);
		w += t;

	}

#ifdef KAT
  {
    int k;
    printf("decrypt e: positions");
    for (k = 0;k < SYS_N;++k)
      if (e[k/8] & (1 << (k&7)))
        printf(" %d",k);
    printf("\n");
  }
#endif

#ifdef TIME_MEASUREMENT
	gettimeofday(&start_synd, NULL);
#endif
#ifdef SYND_KERNEL
	synd_host(s_cmp, g, L, e);
#endif
#ifndef SYND_KERNEL
	synd_sw_host(s_cmp, g, L, e);
#endif
#ifdef TIME_MEASUREMENT
	gettimeofday(&end_synd, NULL);
	get_event_time(&start_synd, &end_synd, &sum_total_synd, &times_total_synd);
#endif

	check = w;
	check ^= SYS_T;

	for (i = 0; i < SYS_T*2; i++)
		check |= s[i] ^ s_cmp[i]; 

	check -= 1;
	check >>= 15;

	return check ^ 1;
}

