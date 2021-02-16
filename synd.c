/*
  This file is for syndrome computation
*/

//#include "synd.h"
//#include "kat_kem.h"
//#include <CL/opencl.h>
//#include <CL/cl_ext.h>
//#include <stdlib.h>
//
//#include "params.h"
//#include "root.h"
//#include <sys/time.h>
//#include <stdio.h>

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

/* input: Goppa polynomial f, support L, received word r */
/* output: out, the syndrome of length 2t */

double sum_synd=0.0;
int times_synd=0;

double sum_synd_2=0.0;
int times_synd_2=0;

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

#ifdef TIME_MEASUREMENT
	cl_event event, event_2;
#endif

	memcpy(ptr_f_in, f, sizeof(gf)*(SYS_T+1));

	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N/2);
	memcpy(ptr_L_in_2, (L+SYS_N/2), sizeof(gf)*SYS_N/2);

	memcpy(ptr_r_in, r, sizeof(unsigned char)*MAT_COLS/2);
	memcpy(ptr_r_in_2, (r+MAT_COLS/2), sizeof(unsigned char)*MAT_COLS/2);

//#ifdef FUNC_CORRECTNESS
//	gf *out_validate = (gf *)malloc(sizeof(gf)*2*SYS_T);
//#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)3, &pt_list_synd[1], 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffers\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)3, &pt_list_synd_2[1], 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffers\n");
    	return EXIT_FAILURE;
    }
	#endif


	#ifdef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_synd, 0, NULL, &event);
    err = clEnqueueTask(commands, kernel_synd_2, 0, NULL, &event_2);
	#endif
	#ifndef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_synd, 0, NULL, NULL);
	#endif
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel synd\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_out_out, CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_out_out_2, CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif


#ifdef TIME_MEASUREMENT
    //probably put n in the place of 1 and wait for all
	clWaitForEvents(1, &event);
	clWaitForEvents(1, &event_2);
#endif
    clFinish(commands);



//    out_xor_red = (gf *)malloc(sizeof(gf)*2*SYS_T);
    for(int i=0; i<2*SYS_T; i++){
    	*(out+i) = *(ptr_out_out+i) ^ *(ptr_out_out_2+i);
    }


//	#ifdef FUNC_CORRECTNESS
////    gf tmp;
//	gf validate_mat[SYS_N];
//	synd_sw_host(out_validate, f, L, r);
//	for (int i=0;i<2*SYS_T;i++){
////		tmp = *(ptr_out_out+i)^*(ptr_out_out_2+i);
//		if (*(out_validate+i) != *(out+i)){\
//			printf("\nERROR in %d: Expected %d, got %d\n", i, *(out_validate+i), *(out+i));
//		}
//	}
//	#endif


	#ifdef TIME_MEASUREMENT
	cl_ulong time_start, time_start_2;
	cl_ulong time_end, time_end_2;

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	clGetEventProfilingInfo(event_2, CL_PROFILING_COMMAND_START, sizeof(time_start_2), &time_start_2, NULL);
	clGetEventProfilingInfo(event_2, CL_PROFILING_COMMAND_END, sizeof(time_end_2), &time_end_2, NULL);

	double nanoSeconds = time_end-time_start;
	sum_synd += nanoSeconds;
	times_synd = times_synd + 1;
//	printf("Syndr kernel: OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);

	double nanoSeconds_2 = time_end_2-time_start_2;
	sum_synd_2 += nanoSeconds_2;
	times_synd_2 = times_synd_2 + 1;
//	printf("Syndr kernel_2: OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds_2 / 1000000.0);
	#endif


}
#endif

