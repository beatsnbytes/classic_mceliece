/*
  This file is for Niederreiter decryption
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
#include "kat_kem.h"
#include "crypto_kem.h"

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

	unsigned char r[ SYS_N/8 ]; //should it be defined only inside the init function?
	gf g[ SYS_T+1 ]; //should it be defined only inside the init function?

	gf L[ SYS_N ];

	gf s[ SYS_T*2 ];
	gf s_cmp[ SYS_T*2 ];
	gf locator[ SYS_T+1 ];
	gf images[ SYS_N ];

//	gf t;//Probably dont need it as t is defined inside the function now



#ifdef INIT_KERNEL
	init_host(g, r, c, sk);
#endif
#ifndef INIT_KERNEL
	init_sw_host(g, r, c, sk);
#endif


#ifdef SUPPORT_KERNEL
	support_gen_host(L, sk);
#endif
#ifndef SUPPORT_KERNEL
	support_gen_sw_host(L, sk);
#endif

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

#ifdef BM_KERNEL
	bm_host(locator, s);
#endif
#ifndef BM_KERNEL
	bm_sw_host(locator, s);
#endif

#ifdef ROOT_KERNEL
	root_host(images, locator, L);
#endif
#ifndef ROOT_KERNEL
	root_sw_host(images, locator, L);
#endif

	
//	for (i = 0; i < SYS_N/8; i++)
//		e[i] = 0;
//
//	for (i = 0; i < SYS_N; i++)
//	{
//		t = gf_iszero(images[i]) & 1;
//
//		e[ i/8 ] |= t << (i%8);
//		w += t;
//	}

	//TODO remove & if I finally pass the output w as scalar
	#ifdef MID_COMPUTE
		mid_host(&w, e, images);
	#endif
	#ifndef MID_COMPUTE
		mid_sw_host(&w, e, images);
	#endif


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

//	check = w;
//	check ^= SYS_T;
//
//	for (i = 0; i < SYS_T*2; i++)
//		check |= s[i] ^ s_cmp[i];
//
//	check -= 1;
//	check >>= 15;
//
//	return check ^ 1;


#ifdef CHECK
	check_host(&check, &w, s, s_cmp);
#endif
#ifndef CHECK
	check_sw_host(&check, &w, s, s_cmp);
#endif

	//The ^1 has happened inside the kernel
	return check;


}

#ifdef INIT_KERNEL

void init_host(gf *g, unsigned char *r, const unsigned char *c, unsigned char *sk)
{


	cl_event events_enq[2], event_migr_tohost[2], event_migr_tokern[1];


	memcpy(ptr_cinit_in, c, sizeof(unsigned char)*(SYND_BYTES));
	memcpy(ptr_skinit_in, sk, sizeof(unsigned char)*(crypto_kem_SECRETKEYBYTES));

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_cinit_in, 0, 0, NULL, &event_migr_tokern[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer cinit\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_skinit_in, 0, 0, NULL, &event_migr_tokern[1]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer cinit\n");
    	return EXIT_FAILURE;
    }
	#endif


//	#ifdef TIME_MEASUREMENT
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

	err = clEnqueueTask(commands, init_kernels_list[0], 2, &event_migr_tokern, &events_enq[0]);

//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(init_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_init_kernels, &times_init_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel init\n");
    	return EXIT_FAILURE;
    }
	#endif

    //can I migrate both output matrices with one call? probably yes
    //try with 2 initially to make sure and then make a list of 2 buffers
	err = clEnqueueMigrateMemObjects(commands, init_kernels, &buffer_ginit_out, CL_MIGRATE_MEM_OBJECT_HOST, init_kernels, &events_enq[0], &event_migr_tohost[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_ginit_out\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, init_kernels, &buffer_rinit_out, CL_MIGRATE_MEM_OBJECT_HOST, init_kernels, &events_enq[0], &event_migr_tohost[1]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_rinit_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    //should that be &event_migr_tohost[0]?
    clWaitForEvents(2, &event_migr_tohost);

    memcpy(g, ptr_ginit_out, sizeof(gf)*(SYS_T+1));
    memcpy(r, ptr_rinit_out, sizeof(gf)*(SYS_N/8));



#ifdef FUNC_CORRECTNESS

#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_init_tokern, &times_init_tokern);
//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_init_kernel, &times_init);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_init_tohost, &times_init_tohost);
//#endif

}

#endif

void init_sw_host(gf *g, unsigned char *r, const unsigned char *c, unsigned char *sk)
{

	//probbly not needed to be defined also here
//	unsigned char r[ SYS_N/8 ];
//	gf g[ SYS_T+1 ];

	for (int i = 0; i < SYND_BYTES; i++)       r[i] = c[i];
	for (int i = SYND_BYTES; i < SYS_N/8; i++) r[i] = 0;

	for (int i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; } g[ SYS_T ] = 1;

}


#ifdef MID_KERNEL

void mid_host(int *w, unsigned char *e, gf *images)
{


	cl_event events_enq[1], event_migr_tohost[2], event_migr_tokern[1];


	memcpy(ptr_imagesmid_in, images, sizeof(gf)*(SYS_N));

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_imagesmid_in, 0, 0, NULL, &event_migr_tokern[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer imagesmid\n");
    	return EXIT_FAILURE;
    }
	#endif


//	#ifdef TIME_MEASUREMENT
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

	err = clEnqueueTask(commands, mid_kernels_list[0], 1, &event_migr_tokern, &events_enq[0]);

//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(init_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_init_kernels, &times_init_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel mid\n");
    	return EXIT_FAILURE;
    }
	#endif

    //can I migrate both output matrices with one call? probably yes
    //try with 2 initially to make sure and then make a list of 2 buffers
	err = clEnqueueMigrateMemObjects(commands, mid_kernels, &buffer_wmid_out, CL_MIGRATE_MEM_OBJECT_HOST, mid_kernels, &events_enq[0], &event_migr_tohost[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_wmid_out\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, mid_kernels, &buffer_emid_out, CL_MIGRATE_MEM_OBJECT_HOST, mid_kernels, &events_enq[0], &event_migr_tohost[1]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_emid_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    //should that be &event_migr_tohost[0]?
    clWaitForEvents(2, &event_migr_tohost);

    memcpy(w, ptr_wmid_out, sizeof(int));
    memcpy(e, ptr_emid_out, sizeof(unsigned char)*(SYS_N/8));



#ifdef FUNC_CORRECTNESS

#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_mid_tokern, &times_mid_tokern);
//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_mid_kernel, &times_mid);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_mid_tohost, &times_mid_tohost);
//#endif

}

#endif


void mid_sw_host(int *w, unsigned char *e, gf *images)
{

	gf t;

	for (int i = 0; i < SYS_N/8; i++)
		e[i] = 0;

	for (int i = 0; i < SYS_N; i++)
	{
		t = gf_iszero(images[i]) & 1;

		e[ i/8 ] |= t << (i%8);
		w += t;
	}

}


#ifdef CHECK_KERNEL

void check_host(uint16_t *check, int *w, gf *s, gf *s_cmp)
{

	cl_event events_enq[1], event_migr_tohost[1], event_migr_tokern[3];


	memcpy(ptr_wcheck_in, w, sizeof(int));
	memcpy(ptr_scheck_in, s, sizeof(gf)*(SYS_T*2));
	memcpy(ptr_s_cmpcheck_in, s_cmp, sizeof(gf)(SYS_T*2));

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_wcheck_in, 0, 0, NULL, &event_migr_tokern[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer wcheck\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_scheck_in, 0, 0, NULL, &event_migr_tokern[1]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer scheck\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_s_cmpcheck_in, 0, 0, NULL, &event_migr_tokern[2]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer s_cmpcheck\n");
    	return EXIT_FAILURE;
    }
	#endif




//	#ifdef TIME_MEASUREMENT
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

	err = clEnqueueTask(commands, check_kernels_list[0], 3, &event_migr_tokern, &events_enq[0]);

//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(init_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_init_kernels, &times_init_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel check\n");
    	return EXIT_FAILURE;
    }
	#endif


	err = clEnqueueMigrateMemObjects(commands, check_kernels, &buffer_check_out, CL_MIGRATE_MEM_OBJECT_HOST, check_kernels, &events_enq[0], &event_migr_tohost[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_check_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    memcpy(check, ptr_check_out, sizeof(uint16_t));



#ifdef FUNC_CORRECTNESS

#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_check_tokern, &times_check_tokern);
//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_check_kernel, &times_check);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_check_tohost, &times_check_tohost);
//#endif




}

#endif

void check_sw_host(uint16_t *check, int *w, gf *s, gf *s_cmp)
{

	*check = w;
	*check ^= SYS_T;

	for (i = 0; i < SYS_T*2; i++)
		*check |= s[i] ^ s_cmp[i];

	*check -= 1;
	*check >>= 15;

	*check^= 1;

}

