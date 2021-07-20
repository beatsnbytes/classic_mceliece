/*
  This file is for the Berlekamp-Massey algorithm
  see http://crypto.stanford.edu/~mironov/cs359/massey.pdf
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#include "params.h"
#include "gf.h"
#include "bm.h"
#include "kat_kem.h"

#define min(a, b) ((a < b) ? a : b)

/* the Berlekamp-Massey algorithm */
/* input: s, sequence of field elements */
/* output: out, minimal polynomial of s */
void bm_sw_host(gf *out, gf *s)
{
	int i;

	uint16_t N = 0;
	uint16_t L = 0;
	uint16_t mle;
	uint16_t mne;

	gf T[ SYS_T+1  ];
	gf C[ SYS_T+1 ];
	gf B[ SYS_T+1 ];

	gf b = 1, d, f;

	//

	for (i = 0; i < SYS_T+1; i++)
		C[i] = B[i] = 0;

	B[1] = C[0] = 1;

	//

	for (N = 0; N < 2 * SYS_T; N++)
	{
		d = 0;

		for (i = 0; i <= min(N, SYS_T); i++)
			d ^= gf_mul(C[i], s[ N-i]);
	
		mne = d; mne -= 1;   mne >>= 15; mne -= 1;
		mle = N; mle -= 2*L; mle >>= 15; mle -= 1;
		mle &= mne;

		for (i = 0; i <= SYS_T; i++)			
			T[i] = C[i];

		f = gf_frac(b, d);

		for (i = 0; i <= SYS_T; i++)			
			C[i] ^= gf_mul(f, B[i]) & mne;

		L = (L & ~mle) | ((N+1-L) & mle);

		for (i = 0; i <= SYS_T; i++)			
			B[i] = (B[i] & ~mle) | (T[i] & mle);

		b = (b & ~mle) | (d & mle);

		for (i = SYS_T; i >= 1; i--) B[i] = B[i-1];
		B[0] = 0;
	}

	for (i = 0; i <= SYS_T; i++)
		out[i] = C[ SYS_T-i ];
}

#ifdef BM_KERNEL
void bm_host(gf *out, gf *s)
{

	cl_event events_enq[3], event_migr_tohost, event_migr_tokern[2];


	memcpy(ptr_sbm_in, s, sizeof(gf)*(SYS_T*2));

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_sbm_in, 0, 0, NULL, &event_migr_tokern[0]);
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

    for (int i=0; i<bm_kernels; i++){
    	err = clEnqueueTask(commands, bm_kernels_list[i], 1, &event_migr_tokern, &events_enq[i]);
    }
//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(root_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_root_kernels, &times_root_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel bm\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, bm_kernels, &buffer_outbm_out, CL_MIGRATE_MEM_OBJECT_HOST, bm_kernels, &events_enq[0], &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_outbm_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    memcpy(out, ptr_outbm_out, sizeof(gf)*(SYS_T+1));



#ifdef FUNC_CORRECTNESS

#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_bm_tokern, &times_bm_tokern);
//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_bm_kernel, &times_bm);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_bm_tohost, &times_bm_tohost);
//#endif

}
#endif

