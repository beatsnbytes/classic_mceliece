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
	cl_event event;
#endif

	memcpy(ptr_f_in, f, sizeof(gf)*(SYS_T+1));
	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N);
	memcpy(ptr_r_in, r, sizeof(unsigned char)*MAT_COLS);

#ifdef FUNC_CORRECTNESS
	gf *out_validate = (gf *)malloc(sizeof(gf)*2*SYS_T);
#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)3, &pt_list_synd[1], 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffers\n");
    	return EXIT_FAILURE;
    }
	#endif

	#ifdef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_synd, 0, NULL, &event);
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

#ifdef TIME_MEASUREMENT
	clWaitForEvents(1, &event);
#endif
    clFinish(commands);




	#ifdef FUNC_CORRECTNESS
    gf validate_mat[SYS_N];
    synd_sw_host(out_validate, f, L, r);
    for (int i=0;i<2*SYS_T;i++){
        if (*(out_validate+i) != *(ptr_r_out+i)){\
        	printf("\nERROR in %d: Expected %d, got %d\n", i, *(out_validate+i), *(out+i));
        }
    }
	#endif

    memcpy(out, ptr_out_out, sizeof(gf)*2*SYS_T);

	#ifdef TIME_MEASUREMENT
	cl_ulong time_start;
	cl_ulong time_end;

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	double nanoSeconds = time_end-time_start;
	sum_synd += nanoSeconds;
	times_synd = times_synd + 1;
	#endif


}
#endif

