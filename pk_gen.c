/*
  This file is for public-key generation
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#include "kat_kem.h"
#include "controlbits.h"
#include "uint64_sort.h"
#include "pk_gen.h"
#include "params.h"
#include "benes.h"
#include "root.h"
#include "util.h"

double sum_elim=0.0;
int times_elim=0;

#ifdef GAUSSIAN_ELIMINATION_KERNEL
void gaussian_elimination_host(unsigned char mat[ GFBITS * SYS_T ][ SYS_N/8 ]) {

	#ifdef TIME_MEASUREMENT
	cl_event event;
	#endif
//	memcpy(ptr_mat_in, mat, sizeof(unsigned char) * MAT_SIZE);
//
//	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_in, 0, 0, NULL, NULL);
//	#ifdef OCL_API_DEBUG
//    if (err != CL_SUCCESS) {
//    	printf("FAILED to enqueue buffer_mat\n");
//    	return EXIT_FAILURE;
//    }
//	#endif

	#ifdef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_gaussian_elimination, 0, NULL, &event);
	#endif
	#ifndef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_gaussian_elimination, 0, NULL, NULL);
	#endif
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_out, CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
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



    //TODO Fix func correctness of gauss elim
	#ifdef FUNC_CORRECTNESS
//    unsigned char *validate_mat = (unsigned char *)malloc(MAT_ROWS * MAT_COLS * sizeof(unsigned char));
    unsigned char validate_mat[MAT_ROWS][MAT_COLS];
	for (int i=0;i<MAT_ROWS;i++){
		for(int j=0;j<MAT_COLS;j++){
			validate_mat[i][j] = *(*mat+i*MAT_COLS+j);
		}
	}

	gaussian_elimination_sw(validate_mat);
	for (int i=0;i<MAT_ROWS;i++){
		for(int j=0;j<MAT_COLS;j++){
			if (validate_mat[i][j] != *(ptr_mat_out+i*MAT_COLS+j)){\
				printf("\nERROR: Expected %d, got %d\n", validate_mat[i][j], *(ptr_mat_out+i*MAT_COLS+j));
			}
		}
	}
	#endif

    memcpy(mat, ptr_mat_out, sizeof(unsigned char) * MAT_SIZE);


	#ifdef TIME_MEASUREMENT
	cl_ulong time_start;
	cl_ulong time_end;

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	double nanoSeconds = time_end-time_start;
	sum_elim += nanoSeconds;
	times_elim = times_elim + 1;
	#endif


}
#endif




// Code that contains the call to gaussian elimination hw kernel
#ifdef GAUSSIAN_ELIMINATION_KERNEL
int pk_gen_host(unsigned char * pk, unsigned char * sk, uint32_t * perm, int16_t * pi)
{
	int i, j, k;
	int row, c;

	uint64_t buf[ 1 << GFBITS ];

	//unsigned char mat[ PK_NROWS ][ SYS_N/8 ];
	unsigned char mask;
	unsigned char b;

	gf g[ SYS_T+1 ]; // Goppa polynomial
	gf L[ SYS_N ]; // support
	gf inv[ SYS_N ];

	//Allocate 4KB aligned memory for mat
//	unsigned char *mat = (unsigned char *)malloc(MAT_ROWS * MAT_COLS * sizeof(unsigned char));

	g[ SYS_T ] = 1;

	for (i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; }

	for (i = 0; i < (1 << GFBITS); i++)
	{
		buf[i] = perm[i];
		buf[i] <<= 31;
		buf[i] |= i;
	}

	uint64_sort(buf, 1 << GFBITS);

	for (i = 1; i < (1 << GFBITS); i++)
		if ((buf[i-1] >> 31) == (buf[i] >> 31))
			return -1;

	for (i = 0; i < (1 << GFBITS); i++) pi[i] = buf[i] & GFMASK;
	for (i = 0; i < SYS_N;         i++) L[i] = bitrev(pi[i]);

	// filling the matrix

	root(inv, g, L);

	for (i = 0; i < SYS_N; i++)
		inv[i] = gf_inv(inv[i]);

	for (i = 0; i < PK_NROWS; i++)
	for (j = 0; j < SYS_N/8; j++)
		*(ptr_mat_in + i*MAT_COLS + j) = 0;

	for (i = 0; i < SYS_T; i++)
	{
		for (j = 0; j < SYS_N; j+=8)
		for (k = 0; k < GFBITS;  k++)
		{
			b  = (inv[j+7] >> k) & 1; b <<= 1;
			b |= (inv[j+6] >> k) & 1; b <<= 1;
			b |= (inv[j+5] >> k) & 1; b <<= 1;
			b |= (inv[j+4] >> k) & 1; b <<= 1;
			b |= (inv[j+3] >> k) & 1; b <<= 1;
			b |= (inv[j+2] >> k) & 1; b <<= 1;
			b |= (inv[j+1] >> k) & 1; b <<= 1;
			b |= (inv[j+0] >> k) & 1;

			*(ptr_mat_in + (i*GFBITS + k)*MAT_COLS +j/8) = b;
		}

		for (j = 0; j < SYS_N; j++)
			inv[j] = gf_mul(inv[j], L[j]);

	}

	//
	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_in, 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue buffer_mat\n");
    	return EXIT_FAILURE;
    }
	#endif
	//

	// gaussian elimination
	gaussian_elimination_host(ptr_mat_in);

	// If the matrix was not systematic repeat the process
	if (*ptr_mat_in == 255){
		return -1;
	}

	for (i = 0; i < PK_NROWS; i++){
		memcpy(pk + i*PK_ROW_BYTES, (ptr_mat_in + i*MAT_COLS) + PK_NROWS/8, PK_ROW_BYTES);
	}

	return 0;
}
#endif

/* input: secret key sk */
/* output: public key pk */
int pk_gen_sw_host(unsigned char * pk, unsigned char * sk, uint32_t * perm, int16_t * pi)
{
	int i, j, k;
	int row, c;

	uint64_t buf[ 1 << GFBITS ];

	unsigned char mat[ PK_NROWS ][ SYS_N/8 ];
	unsigned char mask;
	unsigned char b;

	gf g[ SYS_T+1 ]; // Goppa polynomial
	gf L[ SYS_N ]; // support
	gf inv[ SYS_N ];

	//

	g[ SYS_T ] = 1;

	for (i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; }

	for (i = 0; i < (1 << GFBITS); i++)
	{
		buf[i] = perm[i];
		buf[i] <<= 31;
		buf[i] |= i;
	}

	uint64_sort(buf, 1 << GFBITS);

	for (i = 1; i < (1 << GFBITS); i++)
		if ((buf[i-1] >> 31) == (buf[i] >> 31))
			return -1;

	for (i = 0; i < (1 << GFBITS); i++) pi[i] = buf[i] & GFMASK;
	for (i = 0; i < SYS_N;         i++) L[i] = bitrev(pi[i]);

	// filling the matrix

	root(inv, g, L);

	for (i = 0; i < SYS_N; i++)
		inv[i] = gf_inv(inv[i]);

	for (i = 0; i < PK_NROWS; i++)
	for (j = 0; j < SYS_N/8; j++)
		mat[i][j] = 0;

	for (i = 0; i < SYS_T; i++)
	{
		for (j = 0; j < SYS_N; j+=8)
		for (k = 0; k < GFBITS;  k++)
		{
			b  = (inv[j+7] >> k) & 1; b <<= 1;
			b |= (inv[j+6] >> k) & 1; b <<= 1;
			b |= (inv[j+5] >> k) & 1; b <<= 1;
			b |= (inv[j+4] >> k) & 1; b <<= 1;
			b |= (inv[j+3] >> k) & 1; b <<= 1;
			b |= (inv[j+2] >> k) & 1; b <<= 1;
			b |= (inv[j+1] >> k) & 1; b <<= 1;
			b |= (inv[j+0] >> k) & 1;

			mat[ i*GFBITS + k ][ j/8 ] = b;
		}

		for (j = 0; j < SYS_N; j++)
			inv[j] = gf_mul(inv[j], L[j]);

	}

	// gaussian elimination

	for (i = 0; i < (PK_NROWS + 7) / 8; i++)
	for (j = 0; j < 8; j++)
	{
		row = i*8 + j;

		if (row >= PK_NROWS)
			break;

		for (k = row + 1; k < PK_NROWS; k++)
		{
			mask = mat[ row ][ i ] ^ mat[ k ][ i ];
			mask >>= j;
			mask &= 1;
			mask = -mask;

			for (c = 0; c < SYS_N/8; c++)
				mat[ row ][ c ] ^= mat[ k ][ c ] & mask;
		}

		if ( ((mat[ row ][ i ] >> j) & 1) == 0 ) // return if not systematic
		{
			return -1;
		}

		for (k = 0; k < PK_NROWS; k++)
		{
			if (k != row)
			{
				mask = mat[ k ][ i ] >> j;
				mask &= 1;
				mask = -mask;

				for (c = 0; c < SYS_N/8; c++)
					mat[ k ][ c ] ^= mat[ row ][ c ] & mask;
			}
		}
	}

	for (i = 0; i < PK_NROWS; i++)
		memcpy(pk + i*PK_ROW_BYTES, mat[i] + PK_NROWS/8, PK_ROW_BYTES);

	return 0;
}


