/*
  This file is for Niederreiter encryption
*/

#include "encrypt.h"
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "util.h"
#include "params.h"
#include "randombytes.h"
#include "kat_kem.h"
#include <sys/time.h>
#include "crypto_kem.h"
#include "custom_util.h"

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "gf.h"


double sum_list_syndrome_tokern[1];
double sum_list_syndrome_tohost[1];
double sum_list_syndrome_kernel[8];
int times_syndrome = 0;
int times_syndrome_tohost = 0;
int times_syndrome_tokern = 0;

double sum_syndrome_kernels=0.0;
int times_syndrome_kernels=0;

double sum_total_syndrome=0.0;
int times_total_syndrome=0;

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
static void syndrome_sw_host(unsigned char *s, const unsigned char *pk, unsigned char *e)
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

/* input: public key pk, error vector e */
/* output: syndrome s */
#ifdef SYNDROME_KERNEL
void syndrome_host(unsigned char *s, unsigned char *pk, unsigned char *e)
{

	cl_event events_enq[8], event_migr_tokern, events_migr_tohost;

	memcpy(ptr_pk_in, pk, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES);

	for(int i=0; i<syndrome_kernels; i++){
		memcpy(ptr_e_in_list[i], e, sizeof(unsigned char)*MAT_COLS);
	}


	err = clEnqueueMigrateMemObjects(commands, syndrome_kernels+1, &pt_list_syndrome_combined[0], 0, 0, NULL, &event_migr_tokern);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue pt_list_syndrome\n");
		return EXIT_FAILURE;
	}
	#endif

//	#ifdef TIME_MEASUREMENT
//		clWaitForEvents(1, &event_migr_tokern);
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

	for (int i=0; i<syndrome_kernels; i++){
		err = clEnqueueTask(commands, syndrome_kernels_list[i], 1, &event_migr_tokern, &events_enq[i]);
	}

//	#ifdef TIME_MEASUREMENT
//		clWaitForEvents(syndrome_kernels, &events_enq);
//		gettimeofday(&end_kernel, NULL);
//		get_event_time(&start_kernel, &end_kernel, &sum_syndrome_kernels, &times_syndrome_kernels);
//	#endif

	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to execute kernel\n");
		return EXIT_FAILURE;
	}
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_s_out, CL_MIGRATE_MEM_OBJECT_HOST, syndrome_kernels, &events_enq, &events_migr_tohost);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue bufer_res\n");
		return EXIT_FAILURE;
	}
	#endif

	clWaitForEvents(1, &events_migr_tohost);

	memcpy(s, ptr_s_out, sizeof(unsigned char)*SYND_BYTES);


	#ifdef FUNC_CORRECTNESS
	unsigned char validate_mat[SYND_BYTES];
	syndrome_sw_host(validate_mat, pk, e);
	for (int i=0;i<SYND_BYTES;i++){
		if (validate_mat[i] != *(s+i)){\
			printf("\nERROR in %d: Expected %d, got %d\n", i, validate_mat[i], *(s+i));
		}
	}
	#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_syndrome_tokern, &times_syndrome_tokern);
//	cl_profile_print(&events_enq[0], syndrome_kernels, sum_list_syndrome_kernel, &times_syndrome);
//	cl_profile_print(&events_migr_tohost, 1, sum_list_syndrome_tohost, &times_syndrome_tohost);
//#endif


}
#endif




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

#ifdef TIME_MEASUREMENT
  	struct timeval start_syndrome, end_syndrome;
  	gettimeofday(&start_syndrome, NULL);
#endif

	#ifdef SYNDROME_KERNEL
	syndrome_host(s, pk, e);
	#endif
	#ifndef SYNDROME_KERNEL
	syndrome_sw_host(s, pk, e);
	#endif

#ifdef TIME_MEASUREMENT
	gettimeofday(&end_syndrome, NULL);
	get_event_time(&start_syndrome, &end_syndrome, &sum_total_syndrome, &times_total_syndrome);
#endif

}

