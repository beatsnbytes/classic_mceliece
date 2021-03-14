/*
  This file is for syndrome computation
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "synd.h"
#include "params.h"
#include "root.h"
#include "gf.h"
#include "kat_kem.h"
#include "custom_util.h"

/* input: Goppa polynomial f, support L, received word r */
/* output: out, the syndrome of length 2t */

double sum_list_synd_tokern[1];
double sum_list_synd_tohost[1];
double sum_list_synd_kernel[8];
int times_synd = 0;
int times_synd_tohost = 0;
int times_synd_tokern = 0;

double sum_synd_kernels=0.0;
int times_synd_kernels=0;




void synd_sw_host(gf *out, gf* f , gf *L, unsigned char *r)
{
	int i, j;
	gf e, e_inv, c;

	for (j = 0; j < 2*SYS_T; j++)
		out[j] = 0;

	for (i = 0; i < SYS_N; i++)
	{
		c = (r[i/8] >> (i%8)) & 1;

		e = eval(f, L[i]);
		e_inv = gf_inv(gf_mul(e,e));

		for (j = 0; j < 2*SYS_T; j++)
		{
			out[j] = gf_add(out[j], gf_mul(e_inv, c));
			e_inv = gf_mul(e_inv, L[i]);
		}
	}
}


#ifdef SYND_KERNEL
void synd_host(gf *out, gf *f, gf *L, unsigned char *r)
{

	cl_event events_enq[8], event_migr_tohost, event_migr_tokern;
	gf xor_reduction;


	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N);
	memcpy(ptr_r_in, r, sizeof(unsigned char)*MAT_COLS);

	for (int i=0; i<synd_kernels; i++){
		memcpy(ptr_f_in_list[i], f, sizeof(gf)*(SYS_T+1));
	}

	err = clEnqueueMigrateMemObjects(commands, 3+(synd_kernels-1), &pt_list_synd_combined[0], 0, 0, NULL, &event_migr_tokern);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffers\n");
    	return EXIT_FAILURE;
    }
	#endif

	#ifdef TIME_MEASUREMENT
		struct timeval start_kernel, end_kernel;
		gettimeofday(&start_kernel, NULL);
	#endif

    for (int i=0; i<synd_kernels; i++){
    	err = clEnqueueTask(commands, synd_kernels_list[i], 1, &event_migr_tokern, &events_enq[i]);
    }
    	#ifdef TIME_MEASUREMENT
    		clWaitForEvents(synd_kernels, &events_enq);
    		gettimeofday(&end_kernel, NULL);
    		get_event_time(&start_kernel, &end_kernel, &sum_synd_kernels, &times_synd_kernels);
    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel synd\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, synd_kernels, &pt_list_synd_combined_out[0], CL_MIGRATE_MEM_OBJECT_HOST, synd_kernels, &events_enq[0], &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    for(int i=0; i<2*SYS_T; i++){
    	xor_reduction = 0;
    	for(int j=0; j<synd_kernels; j++){
    		xor_reduction ^= *(ptr_out_out_list[j] + i);
    	}
       *(out+i) = xor_reduction;
    }




#ifdef FUNC_CORRECTNESS
	gf *out_validate = (gf *)malloc(sizeof(gf)*2*SYS_T);
	synd_sw_host(out_validate, f, L, r);
	for (int i=0;i<2*SYS_T;i++){
		if (*(out_validate+i) != *(out+i)){\
			printf("\nERROR in %d: Expected %d, got %d\n", i, *(out_validate+i), *(out+i));
		}
	}
#endif

#ifdef TIME_MEASUREMENT
	cl_profile_print(&event_migr_tokern, 1, sum_list_synd_tokern, &times_synd_tokern);
	cl_profile_print(&events_enq[0], synd_kernels, sum_list_synd_kernel, &times_synd);
	cl_profile_print(&event_migr_tohost, 1, sum_list_synd_tohost, &times_synd_tohost);
#endif

}
#endif


