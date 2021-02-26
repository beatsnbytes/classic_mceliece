#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gf.h"

void cl_profile_print(cl_event *event, int event_num, double *sum_list, int *times){

	cl_ulong time_queue, time_submit, time_start, time_end;
	double miliSeconds_q_sub, miliSeconds_sub_start, miliSeconds_start_end, miliSeconds_total;

	*times += 1;
	for(int i=0; i<event_num; i++){

		clGetEventProfilingInfo(*(event+i), CL_PROFILING_COMMAND_QUEUED, sizeof(time_queue), &time_queue, NULL);
		clGetEventProfilingInfo(*(event+i), CL_PROFILING_COMMAND_SUBMIT, sizeof(time_submit), &time_submit, NULL);
		clGetEventProfilingInfo(*(event+i), CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
		clGetEventProfilingInfo(*(event+i), CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

		miliSeconds_q_sub = (time_submit - time_queue)/1000000.0;
		miliSeconds_sub_start = (time_start - time_submit)/1000000.0;
		miliSeconds_start_end = (time_end - time_start)/1000000.0;
		miliSeconds_total = (time_end - time_queue)/1000000.0;
		sum_list[i] += miliSeconds_total;

#ifdef ITERATION_PRINT
		printf("Event %d: Queued --> Submit time is: %0.3f milliseconds \n", i, miliSeconds_q_sub);
		printf("Event %d: Submit --> Start time is: %0.3f milliseconds \n", i, miliSeconds_sub_start);
		printf("Event %d: Start --> End time is: %0.3f milliseconds \n", i, miliSeconds_start_end);
		printf("Event %d: Total time is: %0.3f milliseconds \n", i, miliSeconds_total);
#endif
	}

	return;
}

void get_event_time(struct timeval *start, struct timeval *end, double *event_sum, int *event_times){

	long seconds, microseconds;
	double elapsed;

	seconds = end->tv_sec - start->tv_sec;
	microseconds = end->tv_usec - start->tv_usec;
	elapsed = seconds + microseconds*0.000001;
	*event_sum += elapsed*1000.0;
	*event_times += 1;
#ifdef ITERATION_PRINT
	printf("Time elapsed is %0.3f milliseconds\n", elapsed);
#endif


	return;

}

void print_kernel_execution_time(double *sum_list, int *times, int instances){

	for(int i=0; i<instances; i++){
    		printf("Event instance %d : Avg Execution time is: %0.3f miliseconds \n", i, sum_list[i]/(*times));
    }
	return;
}

void print_event_execution_time(double *sum_list, int *times){

	printf("Avg Execution time is: %0.3f miliseconds \n", *sum_list/(*times));

	return;
}
