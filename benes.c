#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#include "util.h"
#include "transpose.h"
#include "params.h"
#include "benes.h"
#include "kat_kem.h"
#include "crypto_kem.h"

/*
  This file is for Benes network related functions

  For the implementation strategy, see
  https://eprint.iacr.org/2017/793.pdf
*/
/* one layer of the benes network */
static void layer(uint64_t * data, uint64_t * bits, int lgs)
{
	int i, j, s;

	uint64_t d;

	s = 1 << lgs;

	for (i = 0; i < 64; i += s*2)
	for (j = i; j < i+s; j++)
	{

		d = (data[j+0] ^ data[j+s]);
		d &= (*bits++);
		data[j+0] ^= d;
		data[j+s] ^= d;
	}
}

/* input: r, sequence of bits to be permuted */
/*        bits, condition bits of the Benes network */
/*        rev, 0 for normal application; !0 for inverse */
/* output: r, permuted bits */
void apply_benes(unsigned char * r, const unsigned char * bits, int rev)
{
	int i;

	const unsigned char *cond_ptr; 
	int inc, low;

	uint64_t bs[64];
	uint64_t cond[64];

	//

	for (i = 0; i < 64; i++)
	{
		bs[i] = load8(r + i*8);
	}

	if (rev == 0) 
	{
		inc = 256;
		cond_ptr = bits;
	}
	else
	{
		inc = -256;
		cond_ptr = bits + (2*GFBITS-2)*256;
	}

	//

	transpose_64x64(bs, bs);

	for (low = 0; low <= 5; low++) 
	{ 
		for (i = 0; i < 64; i++) cond[i] = load4(cond_ptr + i*4);
		transpose_64x64(cond, cond);
		layer(bs, cond, low); 
		cond_ptr += inc; 
	}
	
	transpose_64x64(bs, bs);
	
	for (low = 0; low <= 5; low++) 
	{ 
		for (i = 0; i < 32; i++) cond[i] = load8(cond_ptr + i*8);
		layer(bs, cond, low); 
		cond_ptr += inc; 
	}
	for (low = 4; low >= 0; low--) 
	{ 
		for (i = 0; i < 32; i++) cond[i] = load8(cond_ptr + i*8);
		layer(bs, cond, low); 
		cond_ptr += inc; 
	}

	transpose_64x64(bs, bs);
	
	for (low = 5; low >= 0; low--) 
	{ 
		for (i = 0; i < 64; i++) cond[i] = load4(cond_ptr + i*4);
		transpose_64x64(cond, cond);
		layer(bs, cond, low); 
		cond_ptr += inc; 
	}

	transpose_64x64(bs, bs);

	//

	for (i = 0; i < 64; i++)
	{
		store8(r + i*8, bs[i]);
	}
}

/* input: condition bits c */
/* output: support s */
void support_gen_sw_host(gf * s, const unsigned char *c)
{
	gf a;
	int i, j;
	unsigned char L[ GFBITS ][ (1 << GFBITS)/8 ];

	for (i = 0; i < GFBITS; i++)
		for (j = 0; j < (1 << GFBITS)/8; j++)
			L[i][j] = 0;

	for (i = 0; i < (1 << GFBITS); i++)
	{
		a = bitrev((gf) i);

		for (j = 0; j < GFBITS; j++)
			L[j][ i/8 ] |= ((a >> j) & 1) << (i%8);
	}
			
	for (j = 0; j < GFBITS; j++)
		apply_benes(L[j], c, 0);

	for (i = 0; i < SYS_N; i++)
	{
		s[i] = 0;
		for (j = GFBITS-1; j >= 0; j--)
		{
			s[i] <<= 1;
			s[i] |= (L[j][i/8] >> (i%8)) & 1;
		}
	}
}

void support_gen_host(gf * s, const unsigned char *c)
{


	cl_event events_enq[11], event_migr_tohost, event_migr_tokern[1];


	memcpy(ptr_csupp_in, c, sizeof(unsigned char)*crypto_kem_SECRETKEYBYTES);


	//TODO move to one list and enqueue together
	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_csupp_in, 0, 0, NULL, &event_migr_tokern[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer Lr\n");
    	return EXIT_FAILURE;
    }
	#endif


//	#ifdef TIME_MEASUREMENT
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

    for (int i=0; i<support_kernels; i++){
    	err = clEnqueueTask(commands, support_kernels_list[i], 1, &event_migr_tokern, &events_enq[i]);
    }
//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(root_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_root_kernels, &times_root_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel support\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, support_kernels, &buffer_ssupp_out, CL_MIGRATE_MEM_OBJECT_HOST, support_kernels, &events_enq[0], &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_ssupp_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    memcpy(s, ptr_ssupp_out, sizeof(gf)*SYS_N);



	#ifdef FUNC_CORRECTNESS

	#endif

	//#ifdef TIME_MEASUREMENT
	//	cl_profile_print(&event_migr_tokern, 1, sum_list_support_tokern, &times_support_tokern);
	//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_root_kernel, &times_support);
	//	cl_profile_print(&event_migr_tohost, 1, sum_list_support_tohost, &times_support_tohost);
	//#endif


}
