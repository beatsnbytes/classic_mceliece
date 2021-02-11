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

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

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
	int i, j, eq, count;

	union 
	{
		uint16_t nums[ SYS_T*2 ];
		unsigned char bytes[ SYS_T*2 * sizeof(uint16_t) ];
	} buf;

	uint16_t ind[ SYS_T ];
	unsigned char mask;	
	unsigned char val[ SYS_T ];	

	while (1)
	{
		randombytes(buf.bytes, sizeof(buf));

		for (i = 0; i < SYS_T*2; i++)
			buf.nums[i] = load_gf(buf.bytes + i*2);

		// moving and counting indices in the correct range

		count = 0;
		for (i = 0; i < SYS_T*2 && count < SYS_T; i++)
			if (buf.nums[i] < SYS_N)
				ind[ count++ ] = buf.nums[i];
		
		if (count < SYS_T) continue;

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

	#ifdef TIME_MEASUREMENT
	cl_event event;
	#endif

	memcpy(ptr_pk_in, pk, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES);
	memcpy(ptr_e_in, e, sizeof(unsigned char)*MAT_COLS);



	err = clEnqueueMigrateMemObjects(commands, (cl_uint)3, pt_list_syndrome, 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue pt_list_syndrome\n");
		return EXIT_FAILURE;
	}
	#endif

	#ifdef TIME_MEASUREMENT
	err = clEnqueueTask(commands, kernel_syndrome, 0, NULL, &event);
	#endif
	#ifndef TIME_MEASUREMENT
	err = clEnqueueTask(commands, kernel_syndrome, 0, NULL, NULL);
	#endif
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to execute kernel\n");
		return EXIT_FAILURE;
	}
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &pt_list_syndrome[2], CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue bufer_res\n");
		return EXIT_FAILURE;
	}
	#endif

	#ifdef TIME_MEASUREMENT
	clWaitForEvents(1, &event);
	#endif
	clFinish(commands);


	memcpy(s, ptr_s_out, sizeof(unsigned char)*SYND_BYTES);

	#ifdef FUNC_CORRECTNESS
	unsigned char validate_mat[SYND_BYTES];
	syndrome_sw_host(validate_mat, pk, e);
	for (int i=0;i<SYND_BYTES;i++){
		if (validate_mat[i] != *(s+i)){\
			printf("\nERROR: Expected %d, got %d\n", validate_mat[i], *(s+i));
		}
	}
	#endif


	#ifdef TIME_MEASUREMENT
	cl_ulong time_start;
	cl_ulong time_end;

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	double nanoSeconds = time_end-time_start;
	sum_syndrome += nanoSeconds;
	times_syndrome = times_syndrome + 1;
//	printf("Syndrome kernel: OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);
	#endif

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

	#ifdef SYNDROME_KERNEL
	syndrome_host(s, pk, e);
	#endif
	#ifndef SYNDROME_KERNEL
	syndrome_sw_host(s, pk, e);
	#endif

}

