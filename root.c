/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "root.h"
#include "params.h"
#include "gf.h"
#include "kat_kem.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double sum_eval=0;
int times_eval=0;

/* input: polynomial f and field element a */
/* return f(a) */
void eval_sw_host(gf *f, gf *a, gf* out)
{

	int i, j, k;
	gf r, r_tmp;
	uint32_t tmp;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;
	
	for (k = 0; k < SYS_N; k++)
	{

		r = f[ SYS_T ];
		for (i = SYS_T-1; i >= 0; i--)
		{

			t0 = r;
			t1 = a[k];

			tmp = t0 * (t1 & 1);

			for (j = 1; j < GFBITS; j++)
				tmp ^= (t0 * (t1 & (1 << j)));

			t = tmp & 0x7FC000;;
			tmp ^= (t >> 9) ^ (t >> 12);

			t = tmp & 0x3000;
			tmp ^= (t >> 9) ^ (t >> 12);

			r_tmp = tmp & ((1 << GFBITS)-1);

			r = r_tmp ^ f[i];
		}

		out[k]=r;
	}
}


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


#ifdef EVAL_KERNEL
void eval_host(gf *f, gf *a, gf *out)
{

	#ifdef TIME_MEASUREMENT
	cl_event event;
	#endif
	
	memcpy(ptr_f_in, f, sizeof(gf)*(SYS_T+1));
	memcpy(ptr_a_in, a, sizeof(gf)*SYS_N);



	err = clEnqueueMigrateMemObjects(commands, (cl_uint)3, pt_list_eval, 0, 0, NULL, NULL);
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to enqueue buffer_mat\n");
    	return EXIT_FAILURE;
    }
	#endif

	#ifdef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_eval, 0, NULL, &event);
	#endif
	#ifndef TIME_MEASUREMENT
    err = clEnqueueTask(commands, kernel_eval, 0, NULL, NULL);
	#endif
	#ifdef OCL_API_DEBUG
    if (err != CL_SUCCESS) {
    	printf("FAILED to execute kernel\n");
    	return EXIT_FAILURE;
    }
	#endif

	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &pt_list_eval[2], CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
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


    memcpy(out, ptr_r_out, sizeof(gf)*SYS_N);

	#ifdef FUNC_CORRECTNESS
    gf validate_mat[SYS_N];
    eval_sw_host(f, a, validate_mat);
    for (int i=0;i<SYS_N;i++){
        if (validate_mat[i] != *(out+i)){\
        	printf("\nERROR: Expected %d, got %d\n", validate_mat[i], *(out+i));
        }
    }
	#endif


	#ifdef TIME_MEASUREMENT
	cl_ulong time_start;
	cl_ulong time_end;

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	double nanoSeconds = time_end-time_start;
	sum_eval += nanoSeconds;
	times_eval = times_eval + 1;
	#endif

}
#endif


/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root(gf *out, gf *f, gf *L)
{
	int i;

    for (i = 0; i < SYS_N; i++)
            out[i] = eval(f, L[i]);

}

