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
#include "custom_util.h"

/* input: Goppa polynomial f, support L, received word r */
/* output: out, the syndrome of length 2t */

double sum_list_synd_tokern[1];
double sum_list_synd_tohost[1];
double sum_list_synd_kernel[2];
int times_synd = 0;
int times_synd_tohost = 0;
int times_synd_tokern = 0;

double sum_synd_kernels=0.0;
int times_synd_kernels=0;

double sum_list_synd_last_tokern[1];
double sum_list_synd_last_tohost[1];
double sum_list_synd_last_kernel[1];
int times_synd_last = 0;
int times_synd_last_tohost = 0;
int times_synd_last_tokern = 0;


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
	cl_event events_enq[2], event_migr_tohost, event_migr_tokern;
#endif

	memcpy(ptr_f_in, f, sizeof(gf)*(SYS_T+1));
	memcpy(ptr_f_in_2, f, sizeof(gf)*(SYS_T+1));

	memcpy(ptr_f_in_last, f, sizeof(gf)*(SYS_T+1));

	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N);
	memcpy(ptr_r_in, r, sizeof(unsigned char)*MAT_COLS);


//	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N/2);
//	memcpy(ptr_L_in_2, (L+SYS_N/2), sizeof(gf)*SYS_N/2);
//
//	memcpy(ptr_r_in, r, sizeof(unsigned char)*MAT_COLS/2);
//	memcpy(ptr_r_in_2, (r+MAT_COLS/2), sizeof(unsigned char)*MAT_COLS/2);


#ifdef FUNC_CORRECTNESS
	gf *out_validate = (gf *)malloc(sizeof(gf)*2*SYS_T);
#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)5, &pt_list_synd_combined, 0, 0, NULL, &event_migr_tokern);
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
    err = clEnqueueTask(commands, kernel_synd, 1, &event_migr_tokern, &events_enq[0]);
    err = clEnqueueTask(commands, kernel_synd_2, 1, &event_migr_tokern, &events_enq[1]);

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel synd\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)2, &pt_list_synd_combined_out, CL_MIGRATE_MEM_OBJECT_HOST, 2, &events_enq, &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);


    for(int i=0; i<2*SYS_T; i++){
    	*(out+i) = *(ptr_out_out+i) ^ *(ptr_out_out_2+i);
    }


#ifdef FUNC_CORRECTNESS
	gf validate_mat[SYS_N];
	synd_sw_host(out_validate, f, L, r);
	for (int i=0;i<2*SYS_T;i++){
		if (*(out_validate+i) != *(out+i)){\
			printf("\nERROR in %d: Expected %d, got %d\n", i, *(out_validate+i), *(out+i));
		}
	}
#endif

//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_synd_tokern, &times_synd_tokern);
//	cl_profile_print(&events_enq[0], 2, sum_list_synd_kernel, &times_synd);
//	cl_profile_print(&event_migr_tohost, 1, sum_list_synd_tohost, &times_synd_tohost);
//#endif

}
#endif



#ifdef SYND_KERNEL
void synd_host_last(gf *out, gf *f, gf *L, unsigned char *r)
{

#ifdef TIME_MEASUREMENT
	cl_event events_enq[2], event_migr_tohost, event_migr_tokern;
#endif

	memcpy(ptr_r_in_last, r, sizeof(unsigned char)*MAT_COLS);

	//
	cl_mem temp_list[3];
	memcpy(ptr_L_in, L, sizeof(gf)*SYS_N);
	memcpy(ptr_f_in_last, f, sizeof(gf)*(SYS_T+1));
	temp_list[0]= buffer_f_in_last;
	temp_list[1]= buffer_L_in;
	temp_list[2]= buffer_r_in_last;

	//



#ifdef FUNC_CORRECTNESS
	gf *out_validate = (gf *)malloc(sizeof(gf)*2*SYS_T);
#endif

	err = clEnqueueMigrateMemObjects(commands, 3, &temp_list, 0, 0, NULL, &event_migr_tokern);
//	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_r_in_last, 0, 0, NULL, &event_migr_tokern);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue input buffers\n");
    	return EXIT_FAILURE;
    }
	#endif

    err = clEnqueueTask(commands, kernel_synd_last, 1, &event_migr_tokern, &events_enq[0]);

	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel synd\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, 1, &buffer_out_out_last, CL_MIGRATE_MEM_OBJECT_HOST, 1, &events_enq[0], &event_migr_tohost);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue bufer_res\n");
    	return EXIT_FAILURE;
    }
	#endif

    clWaitForEvents(1, &event_migr_tohost);

    memcpy(out, ptr_out_out_last, sizeof(gf)*2*SYS_T);

#ifdef FUNC_CORRECTNESS
	gf validate_mat[SYS_N];
	synd_sw_host(out_validate, f, L, r);
	for (int i=0;i<2*SYS_T;i++){
		if (*(out_validate+i) != *(out+i)){\
			printf("\nERROR in %d: Expected %d, got %d\n", i, *(out_validate+i), *(out+i));
		}
	}
#endif

#ifdef TIME_MEASUREMENT
	cl_profile_print(&event_migr_tokern, 1, sum_list_synd_last_tokern, &times_synd_last_tokern);
	cl_profile_print(&events_enq[0], 1, sum_list_synd_last_kernel, &times_synd_last);
	cl_profile_print(&event_migr_tohost, 1, sum_list_synd_last_tohost, &times_synd_last_tohost);
#endif

}
#endif





