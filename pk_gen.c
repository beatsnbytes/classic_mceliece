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
#include "crypto_hash.h"
#include "uint64_sort.h"
#include "crypto_kem.h"
#include "pk_gen.h"
#include "params.h"
#include "benes.h"
#include "root.h"
#include "util.h"


double sum_list_elim_tokern[1];
double sum_list_elim_tohost[1];
double sum_list_elim_kernel[1];
double sum_elim=0.0;
int times_elim = 0;
int times_elim_tohost = 0;
int times_elim_tokern = 0;

double sum_pk_loop=0.0;
int times_pk_loop=0;

double sum_parallel=0.0;
int times_parallel=0;


#ifdef GAUSSIAN_ELIMINATION_KERNEL


void parallel_sw_part(unsigned char * sk, unsigned char * seed, unsigned char mat[ GFBITS * SYS_T ][ SYS_N/8 ], int16_t * pi){

	int i, j, k;
	unsigned char *rp;
	unsigned char r[ SYS_N/8 + (1 << GFBITS)*sizeof(uint32_t) + SYS_T*2 + 32 ];
	unsigned char b;
	unsigned char *skp;

	gf f[ SYS_T ]; // element in GF(2^mt)
	gf irr[ SYS_T ]; // Goppa polynomial
	uint32_t perm[ 1 << GFBITS ]; // random permutation as 32-bit integers
	gf g[ SYS_T+1 ]; // Goppa polynomial
	gf L[ SYS_N ]; // support
	gf inv[ SYS_N ];
	uint64_t buf[ 1 << GFBITS ];

	//***** END PART OF KEM_KEYPAIR *****//
	while(1){
		rp = &r[ sizeof(r)-32 ];
		skp = sk;

		// expanding and updating the seed
		shake(r, sizeof(r), seed, 33);

		memcpy(skp, seed+1, 32);
		skp += 32 + 8;
		memcpy(seed+1, &r[ sizeof(r)-32 ], 32);

		// generating irreducible polynomial
		rp -= sizeof(f);

		for (i = 0; i < SYS_T; i++) {
			f[i] = load_gf(rp + i*2);

		}

		if (genpoly_gen(irr, f)){
			continue;
		}

		for (i = 0; i < SYS_T; i++)
			store_gf(skp + i*2, irr[i]);

		skp += IRR_BYTES;

		// generating permutation
		rp -= sizeof(perm);

		for (i = 0; i < (1 << GFBITS); i++){
			perm[i] = load4(rp + i*4);
		}

		//***** END PART OF KEM_KEYPAIR *****//


		//***** PART OF PK_GEN *****//
		sk = skp - IRR_BYTES;

		g[ SYS_T ] = 1;

		for (i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; }

		for (i = 0; i < (1 << GFBITS); i++)
		{
			buf[i] = perm[i];
			buf[i] <<= 31;
			buf[i] |= i;
		}

		uint64_sort(buf, 1 << GFBITS);

		for (i = 1; i < (1 << GFBITS); i++){
			if ((buf[i-1] >> 31) == (buf[i] >> 31)){
				// return -1;
				continue; //No need to exit. Just start from the start of the while(1) loop
			}
		}

		break;
	}

	for (i = 0; i < (1 << GFBITS); i++){ pi[i] = buf[i] & GFMASK;}
	for (i = 0; i < SYS_N;         i++){ L[i] = bitrev(pi[i]);}

	// filling the matrix

	root(inv, g, L);

	for (i = 0; i < SYS_N; i++)
		inv[i] = gf_inv(inv[i]);

	for (i = 0; i < PK_NROWS; i++){
		for (j = 0; j < SYS_N/8; j++){
			mat[i][j] = 0;
		}
	}

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
	//***** END PART OF PK_GEN *****//

}


int gaussian_elim(unsigned char mat[ GFBITS * SYS_T ][ SYS_N/8 ], uint *fail){

	int i, j, k;
	int row, c;
	unsigned char mask;

	*fail=0;
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
			mat[0][0]=255;
			*fail=1;
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

	return 0;

}





void gaussian_elimination_host(unsigned char mat[ GFBITS * SYS_T ][ SYS_N/8 ], uint *fail){
//void gaussian_elimination_host(unsigned char *mat, uint *fail) {

	#ifdef TIME_MEASUREMENT
	cl_event event_enq, event_mig_tokern, event_mig_tohost, event_mig_tohost_buffer;
	#endif


	memcpy(ptr_mat_in, mat, sizeof(unsigned char) * MAT_SIZE);

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_in, 0, 0, NULL, &event_mig_tokern);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue buffer_mat\n");
    	return EXIT_FAILURE;
    }
	#endif

	#ifdef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_gaussian_elimination, 1, &event_mig_tokern, &event_enq);
	#endif
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_fail, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_enq, &event_mig_tohost_buffer);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue buffer success info\n");
		return EXIT_FAILURE;
	}
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_out, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_enq, &event_mig_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_mig_tohost_buffer);
    clWaitForEvents(1, &event_mig_tohost);

    *fail = *ptr_fail;
    memcpy(mat, ptr_mat_out, sizeof(unsigned char) * (MAT_ROWS*PK_ROW_BYTES));

#ifdef TIME_MEASUREMENT
	cl_profile_print(&event_mig_tokern, 1, sum_list_elim_tokern, &times_elim_tokern);
    cl_profile_print(&event_enq, 1, sum_list_elim_kernel, &times_elim);
	cl_profile_print(&event_mig_tohost, 1, sum_list_elim_tohost, &times_elim_tohost);
#endif
}
#endif




// Code that contains the call to gaussian elimination hw kernel
#ifdef GAUSSIAN_ELIMINATION_KERNEL
int pk_gen_host(unsigned char * pk, unsigned char * sk, uint32_t * perm, int16_t * pi, unsigned char * sk_initial,  unsigned char * seed_initial)
{

	int i, j, k;
	uint64_t buf[ 1 << GFBITS ];
	unsigned char b;
	gf g[ SYS_T+1 ]; // Goppa polynomial
	gf L[ SYS_N ]; // support
	gf inv[ SYS_N ];

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
		if ((buf[i-1] >> 31) == (buf[i] >> 31)){
			return -1;
		}

	for (i = 0; i < (1 << GFBITS); i++) pi[i] = buf[i] & GFMASK;
	for (i = 0; i < SYS_N;         i++) L[i] = bitrev(pi[i]);

	// filling the matrix

	root(inv, g, L);

	for (i = 0; i < SYS_N; i++)
		inv[i] = gf_inv(inv[i]);

	for (i = 0; i < PK_NROWS; i++)
	for (j = 0; j < SYS_N/8; j++){
		*(ptr_mat_in + i*MAT_COLS + j) = 0;

	}

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



#ifdef TIME_MEASUREMENT
  	struct timeval start_pk_loop, end_pk_loop;
  	gettimeofday(&start_pk_loop, NULL);
#endif

	uint fail=0;
	gaussian_elimination_host(ptr_mat_in, &fail);
	while(fail==1)
	{
	#ifdef TIME_MEASUREMENT
		struct timeval start_parallel, end_parallel;
		gettimeofday(&start_parallel, NULL);
	#endif

		parallel_sw_part(sk_initial, seed_initial, ptr_mat_in, pi);

	#ifdef TIME_MEASUREMENT
			gettimeofday(&end_parallel, NULL);
			get_event_time(&start_parallel, &end_parallel, &sum_parallel, &times_parallel);
	#endif
		gaussian_elimination_host(ptr_mat_in, &fail);
	}

#ifdef TIME_MEASUREMENT
	    gettimeofday(&end_pk_loop, NULL);
	    get_event_time(&start_pk_loop, &end_pk_loop, &sum_pk_loop, &times_pk_loop);
#endif


//
//	unsigned char * sk_parallel = (unsigned char *)malloc(sizeof(unsigned char) * crypto_kem_SECRETKEYBYTES);;
//	int16_t * pi_parallel= (int16_t *)malloc(sizeof(uint16_t) * (1<<GFBITS));
//
//	memcpy(sk_parallel, sk_initial, sizeof(unsigned char) * crypto_kem_SECRETKEYBYTES);
//	memcpy(pi_parallel, pi, sizeof(uint16_t) * (1<<GFBITS));
//	cl_event event_enq, event_mig_tohost_buffer, event_mig_tokern, event_mig_tohost_mat;
//
//#ifdef TIME_MEASUREMENT
//  	struct timeval start_pk_loop, end_pk_loop;
//  	gettimeofday(&start_pk_loop, NULL);
//#endif
//
//	do
//	{
//
//		//Send data and enqueue kernel
//		err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_in, 0, 0, NULL, &event_mig_tokern);
//		#ifdef OCL_API_DEBUG
//	    if (err != CL_SUCCESS) {
//	    	printf("FAILED to enqueue buffer_mat\n");
//	    	return EXIT_FAILURE;
//	    }
//		#endif
//
//
//	    err = clEnqueueTask(commands, kernel_gaussian_elimination, 1, &event_mig_tokern, &event_enq);
//
//		err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_fail, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_enq, &event_mig_tohost_buffer);
//		#ifdef OCL_API_DEBUG
//	    if (err != CL_SUCCESS) {
//	    	printf("FAILED to enqueue buffer success info\n");
//	    	return EXIT_FAILURE;
//	    }
//		#endif
//
//#ifdef TIME_MEASUREMENT
//  	struct timeval start_parallel, end_parallel;
//  	gettimeofday(&start_parallel, NULL);
//#endif
//		memcpy(sk_initial, sk_parallel, sizeof(unsigned char) * crypto_kem_SECRETKEYBYTES);
//		memcpy(pi, pi_parallel, sizeof(uint16_t) * (1<<GFBITS));
//
//	    parallel_sw_part(sk_parallel, seed_initial, ptr_mat_in, pi_parallel);
//
//#ifdef TIME_MEASUREMENT
//	    gettimeofday(&end_parallel, NULL);
//	    get_event_time(&start_parallel, &end_parallel, &sum_parallel, &times_parallel);
//#endif
//	    clWaitForEvents(1, &event_mig_tohost_buffer);
//
//
//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_mig_tokern, 1, sum_list_elim_tokern, &times_elim_tokern);
//	cl_profile_print(&event_enq, 1, sum_list_elim_kernel, &times_elim);
//	cl_profile_print(&event_mig_tohost_buffer, 1, sum_list_elim_tohost, &times_elim_tohost);
////	if(*ptr_fail!=0){
////		times_elim = times_elim - 1;
////		times_elim_tokern = times_elim_tokern - 1;
////		times_elim_tohost = times_elim_tohost - 1;
////	}
//#endif
//
//	}while(*ptr_fail!=0);
//
//
//#ifdef TIME_MEASUREMENT
//	    gettimeofday(&end_pk_loop, NULL);
//	    get_event_time(&start_pk_loop, &end_pk_loop, &sum_pk_loop, &times_pk_loop);
//#endif
//
//	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_mat_out, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_mig_tohost_buffer, &event_mig_tohost_mat);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to enqueue buffer success info\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	clWaitForEvents(1, &event_mig_tohost_mat);

	for (i = 0; i < PK_NROWS; i++)
		memcpy(pk + i*PK_ROW_BYTES, (ptr_mat_out+i*MAT_COLS) + PK_NROWS/8, PK_ROW_BYTES);


//	free(sk_parallel);
//	free(pi_parallel);

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

