/*
  This file is for evaluating a polynomial at one or more field elements
*/



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#include "root.h"
#include "params.h"
#include "gf.h"
#include "kat_kem.h"

/* input: polynomial f and field element a */
/* return f(a) */
gf eval(gf *f, gf a)
{
	int i;
	gf r;
	
	r = f[ SYS_T ];

	for (i = SYS_T-1; i >= 0; i--)
	{
		r = gf_mul(r, a);
		r = gf_add(r, f[i]);
	}

	return r;
}


/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root_sw_host(gf *out, gf *f, gf *L)
{
	int i, j;

	gf r;

	for (j = 0; j < SYS_N; j++){

		r = f[ SYS_T ];

		for (i = SYS_T-1; i >= 0; i--)
		{
			r = gf_mul(r, L[j]);
			r = gf_add(r, f[i]);
		}

		out[j] = r;
	}


}

#ifdef ROOT_KERNEL
void root_host(gf *out, gf *f, gf *L)
{

	cl_event events_enq[11], event_migr_tohost, event_migr_tokern[2];


	memcpy(ptr_Lr_in, L, sizeof(gf)*SYS_N);
	memcpy(ptr_fr_in, f, sizeof(gf)*(SYS_T+1));


	//TODO move to one list and enqueue together
	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_Lr_in, 0, 0, NULL, &event_migr_tokern[0]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer Lr\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_fr_in, 0, 0, NULL, &event_migr_tokern[1]);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffer fr\n");
    	return EXIT_FAILURE;
    }
	#endif

//	#ifdef TIME_MEASUREMENT
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

    for (int i=0; i<root_kernels; i++){
    	err = clEnqueueTask(commands, root_kernels_list[i], 2, &event_migr_tokern, &events_enq[i]);
    }
//    	#ifdef TIME_MEASUREMENT
//    		clWaitForEvents(root_kernels, &events_enq);
//    		gettimeofday(&end_kernel, NULL);
//    		get_event_time(&start_kernel, &end_kernel, &sum_root_kernels, &times_root_kernels);
//    	#endif

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel root\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, root_kernels, &buffer_out_out, CL_MIGRATE_MEM_OBJECT_HOST, root_kernels, &events_enq[0], &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_out_out\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    memcpy(out, ptr_out_out, sizeof(gf)*SYS_N);



#ifdef FUNC_CORRECTNESS

#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_root_tokern, &times_root_tokern);
//	cl_profile_print(&events_enq[0], synd_kernels, sum_list_root_kernel, &times_root);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_root_tohost, &times_root_tohost);
//#endif

}
#endif



