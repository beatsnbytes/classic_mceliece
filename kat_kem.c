/*
   PQCgenKAT_kem.c
   Created by Bassham, Lawrence E (Fed) on 8/29/17.
   Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
   + mods from djb: see KATNOTES
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"
#include "crypto_kem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "params.h"
#include <stdbool.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include <math.h>

#include "pk_gen.h"
#include "decrypt.h"
#include "encrypt.h"
#include "root.h"
#include "gf.h"
#include "controlbits.h"
#include "synd.h"
#include "operations.h"
#include "custom_util.h"

#include <sys/time.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include <valgrind/callgrind.h>

#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_CRYPTO_FAILURE  -4

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

void	fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);

unsigned char entropy_input[48];
unsigned char seed[KATNUM][48];

//GLOBAL VARIABLES
cl_int err;
const cl_uint DATA_SIZE = 4096;
cl_platform_id platform_id;
cl_device_id device_id;
cl_context context;
cl_command_queue commands;
cl_program program;

#ifdef GAUSSIAN_ELIMINATION_KERNEL
cl_kernel kernel_gaussian_elimination;

cl_mem buffer_mat_in;
cl_mem buffer_mat_out;
unsigned char *ptr_mat_in;
unsigned char *ptr_mat_out;
unsigned int *ptr_fail;
unsigned char *success_info_host_ptr;
cl_buffer_region region_success_info;
cl_mem buffer_success_info;
cl_mem buffer_fail;

#endif


#ifdef SYNDROME_KERNEL
int syndrome_kernels = 8 ;
cl_kernel syndrome_kernels_list[8];

const char *syndrome_kernels_name_list[15] = {"syndrome_kernel",
										"syndrome_kernel2_1",
										"syndrome_kernel2_2",
										"syndrome_kernel4_1",
										"syndrome_kernel4_2",
										"syndrome_kernel4_3",
										"syndrome_kernel4_4",
										"syndrome_kernel8_1",
										"syndrome_kernel8_2",
										"syndrome_kernel8_3",
										"syndrome_kernel8_4",
										"syndrome_kernel8_5",
										"syndrome_kernel8_6",
										"syndrome_kernel8_7",
										"syndrome_kernel8_8"
										};


cl_mem pt_list_syndrome_combined[9];

cl_mem buffer_pk_in;
unsigned char *ptr_pk_in;

cl_mem buffer_s_out;
unsigned char *ptr_s_out;

cl_mem buffer_e_in_list[8];
unsigned char *ptr_e_in_list[8];

cl_mem pt_list_syndrome_combined[9];

#endif

#ifdef SYND_KERNEL
int synd_kernels = 16;

cl_kernel synd_kernels_list[16];
const char *synd_kernels_name_list[61] = {"synd_kernel1_1",
										"synd_kernel2_1", "synd_kernel2_2",
										"synd_kernel4_1", "synd_kernel4_2", "synd_kernel4_3", "synd_kernel4_4",
										"synd_kernel8_1", "synd_kernel8_2", "synd_kernel8_3", "synd_kernel8_4","synd_kernel8_5", "synd_kernel8_6", "synd_kernel8_7", "synd_kernel8_8",
										"synd_kernel12_1", "synd_kernel12_2", "synd_kernel12_3", "synd_kernel12_4", "synd_kernel12_5", "synd_kernel12_6", "synd_kernel12_7", "synd_kernel12_8", "synd_kernel12_9","synd_kernel12_10","synd_kernel12_11","synd_kernel12_12",
										"synd_kernel16_1", "synd_kernel16_2", "synd_kernel16_3", "synd_kernel16_4", "synd_kernel16_5", "synd_kernel16_6", "synd_kernel16_7", "synd_kernel16_8", "synd_kernel16_9","synd_kernel16_10","synd_kernel16_11","synd_kernel16_12", "synd_kernel16_13","synd_kernel16_14","synd_kernel16_15","synd_kernel16_16",
										"synd_kernel18_1", "synd_kernel18_2", "synd_kernel18_3", "synd_kernel18_4", "synd_kernel18_5", "synd_kernel18_6", "synd_kernel18_7", "synd_kernel18_8", "synd_kernel18_9","synd_kernel18_10","synd_kernel18_11","synd_kernel18_12", "synd_kernel18_13","synd_kernel18_14","synd_kernel18_15","synd_kernel18_16","synd_kernel18_17","synd_kernel18_18"
										};

cl_mem buffer_out_out_list[16];
gf * ptr_out_out_list[16];
cl_mem buffer_f_in_list[16];
gf *ptr_f_in_list[16];

cl_mem buffer_L_in;
gf *ptr_L_in;

cl_mem buffer_r_in;
unsigned char *ptr_r_in;

cl_mem pt_list_synd_combined[18];
cl_mem pt_list_synd_combined_out[16];


#endif


double sum_keygen=0;
double sum_enc=0;
double sum_dec=0;
int times_keygen=0;
int times_enc=0;
int times_dec=0;

char cl_platform_vendor[1001];
const char* target_device_name ="zcu102_base";
cl_platform_id platforms[16];
cl_uint platform_count;
cl_uint platform_found = 0;
cl_uint num_devices;
cl_uint device_found = 0;
cl_device_id devices[16];
char cl_device_name[1001];
cl_int status;

cl_uint load_file_to_memory(const char *filename, char **result)
{
    cl_uint size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) {
        *result = NULL;
        return -1; // -1 means file opening fail
    }
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) {
        free(*result);
        return -2; // -2 means file reading fail
    }
    fclose(f);
    (*result)[size] = 0;
    return size;
}

void	fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);


unsigned char entropy_input[48];
unsigned char seed[KATNUM][48];

int
main(int argc, char* argv[])
{
	//OpenCL Initialization code STARTS
	// ------------------------------------------------------------------------------------
	// Step 1: Get All PLATFORMS, then search for Target_Platform_Vendor (CL_PLATFORM_VENDOR)
	// ------------------------device_------------------------------------------------------------

	// Get the number of platforms
	// ..................................................

    err = clGetPlatformIDs(16, platforms, &platform_count);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("Error: Failed to find an OpenCL platform!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	printf("INFO: Found %d platforms\n", platform_count);
	#endif



	  // ....................................................................................
	  // step 1:  Search for Platform (ex: Xilinx) using: CL_PLATFORM_VENDOR = Target_Platform_Vendor
	  // Check if the current platform matches Target_Platform_Vendor
	  // ................device_....................................................................

	for (cl_uint iplat=0; iplat<platform_count; iplat++) {
		err = clGetPlatformInfo(platforms[iplat], CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor,NULL);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		if (strcmp(cl_platform_vendor, "Xilinx") == 0) {
			printf("INFO: Selected platform %d from %s\n", iplat, cl_platform_vendor);
			platform_id = platforms[iplat];
			platform_found = 1;
		}
		#endif
	}
	#ifdef OCL_API_DEBUG
	if (!platform_found) {
		printf("ERROR: Platform Xilinx not found. Exit.\n");
		return EXIT_FAILURE;
	}
	#endif

	   // ------------------------------------------------------------------------------------
	   // Step 1:  Get All Devices for selected platform Target_Platform_ID
	   //            then search for Xilinx platform (CL_DEVICE_TYPE_ACCELERATOR = Target_Device_Name)
	   // ------------------------------------------------------------------------------------


    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 16, devices, &num_devices);
	#ifdef OCL_API_DEBUG
	printf("INFO: Found %d devices\n", num_devices);
	if (err != CL_SUCCESS) {
		printf("ERROR: Failed to create a device group!\n");
		printf("ERROR: Test failed\n");
		return -1;
	}
	#endif
	 // ------------------------------------------------------------------------------------
	 // Step 1:  Search for CL_DEVICE_NAME = Target_Device_Name
	 // ............................................................................

   for (cl_uint i=0; i<num_devices; i++) {
		err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 1024, cl_device_name, 0);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("Error: Failed to get device name for device %d!\n", i);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("CL_DEVICE_NAME %s\n", cl_device_name);
		#endif



	   // ............................................................................
	   // Step 1: Check if the current device matches Target_Device_Name
	   // ............................................................................

	   if(strcmp(cl_device_name, target_device_name) == 0) {
			device_id = devices[i];
			device_found = 1;
			#ifdef OCL_API_DEBUG
			printf("Selected %s as the target device\n", cl_device_name);
			#endif
	   }
	}


	// ------------------------------------------------------------------------------------
	// Step 1: Create Context
	// ------------------------------------------------------------------------------------
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (!context) {
		printf("Error: Failed to create a compute context!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	// ------------------------------------------------------------------------------------
	// Step 1: Create Command Queue
	// ------------------------------------------------------------------------------------
	commands = clCreateCommandQueue(context, device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err);
	#ifdef OCL_API_DEBUG
	if (!commands) {
		printf("Error: Failed to create a command commands!\n");
		printf("Error: code %i\n",err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif



   cl_int status;
   unsigned char *kernelbinary;
   char *xclbin = argv[1];

   // ------------------------------------------------------------------
   // Step 1: Load Binary File from a disk to Memory
   // ------------------------------------------------------------------
	#ifdef OCL_API_DEBUG
    printf("INFO: loading xclbin %s\n", xclbin);
	#endif
    cl_uint n_i0 = load_file_to_memory(xclbin, (char **) &kernelbinary);
	#ifdef OCL_API_DEBUG
    if (n_i0 < 0) {
	    printf("failed to load kernel from xclbin: %s\n", xclbin);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
    }
	#endif



	// ------------------------------------------------------------
	// Step 1: Create a program using a Binary File
	// ------------------------------------------------------------
	size_t n0 = n_i0;
	program = clCreateProgramWithBinary(context, 1, &device_id, &n0, (const unsigned char **) &kernelbinary, &status, &err);
	free(kernelbinary);

	// ============================================================================
	// Step 2: Create Program and Kernels
	// ============================================================================
	//   o) Build a Program from a Binary File
	//   o) Create Kernels
	// ============================================================================
	#ifdef OCL_API_DEBUG
	if ((!program) || (err!=CL_SUCCESS)) {
		printf("Error: Failed to create compute program from binary %d!\n", err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	// -------------------------------------------------------------
	// Step 2: Build (compiles and links) a program executable from binary
	// -------------------------------------------------------------
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];

		printf("Error: Failed to build program executable!\n");
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif
	// -------------------------------------------------------------
	// Step 2: Create Kernels
	// -------------------------------------------------------------


#ifdef GAUSSIAN_ELIMINATION_KERNEL
	kernel_gaussian_elimination = clCreateKernel(program, "gaussian_elimination_kernel", &err);
	#ifdef OCL_API_DEBUG
	if (!kernel_gaussian_elimination || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel_gaussian_elimination!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif


	buffer_mat_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char) * MAT_SIZE, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	buffer_mat_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char) * MAT_SIZE, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	buffer_fail = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned int), NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_gaussian_elimination, 0, sizeof(cl_mem), &buffer_mat_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_gaussian_elimination, 1, sizeof(cl_mem), &buffer_mat_out);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_gaussian_elimination, 2, sizeof(cl_mem), &buffer_fail);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_fail");
		return EXIT_FAILURE;
	}
	#endif


	ptr_mat_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_mat_in, false, CL_MAP_WRITE, 0, sizeof(unsigned char) * MAT_SIZE, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	ptr_mat_out = (unsigned char *) clEnqueueMapBuffer(commands, buffer_mat_out, true, CL_MAP_READ, 0, sizeof(unsigned char) * MAT_SIZE, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	ptr_fail = (unsigned int *) clEnqueueMapBuffer(commands, buffer_fail, false, CL_MAP_READ, 0, sizeof(unsigned int), 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_a");
		return EXIT_FAILURE;
	}
	#endif



#endif


#ifdef SYNDROME_KERNEL


		buffer_pk_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to create buffer_pk_in");
			return EXIT_FAILURE;
		}
		#endif


		buffer_s_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*SYND_BYTES, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to create buffer_s_out");
			return EXIT_FAILURE;
		}
		#endif


		ptr_pk_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_pk_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES, 0, NULL, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("ERROR : %d\n", err);
			printf("FAILED to enqueue map buffer_pk_in");
			return EXIT_FAILURE;
		}
		#endif


		ptr_s_out = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out, false, CL_MAP_READ, 0, sizeof(unsigned char)*SYND_BYTES, 0, NULL, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("ERROR : %d\n", err);
			printf("FAILED to enqueue map buffer_s_out");
			return EXIT_FAILURE;
		}
		#endif

		//Iterative generation of buffers/pointer for the rest of the kernels defined
		int syndrome_idx;
		switch(syndrome_kernels){
			case 8:
				syndrome_idx = 7;
				break;
			default:
				break;
		}



		for(int i=0; i<syndrome_kernels; i++){

			int index = syndrome_idx+i;
			syndrome_kernels_list[i] = clCreateKernel(program, syndrome_kernels_name_list[index], &err);
			#ifdef OCL_API_DEBUG
			if (!syndrome_kernels_list[i] || err != CL_SUCCESS) {
				printf("Error: Failed to create compute kernel_syndrome!\n");
				printf("Test failed\n");
				return EXIT_FAILURE;
			}
			#endif


			buffer_e_in_list[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to create buffer_e_in");
				return EXIT_FAILURE;
			}
			#endif

			err = clSetKernelArg(syndrome_kernels_list[i], 0, sizeof(cl_mem), &buffer_pk_in);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_pk_in");
				return EXIT_FAILURE;
			}
			#endif

			err = clSetKernelArg(syndrome_kernels_list[i], 1, sizeof(cl_mem), &buffer_e_in_list[i]);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_e_in");
				return EXIT_FAILURE;
			}
			#endif


			err = clSetKernelArg(syndrome_kernels_list[i], 2, sizeof(cl_mem), &buffer_s_out);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_s_out");
				return EXIT_FAILURE;
			}
			#endif


			ptr_e_in_list[i] = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_in_list[i], false, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("ERROR : %d\n", err);
				printf("FAILED to enqueue map buffer_e_in");
				return EXIT_FAILURE;
			}
			#endif

		}


		pt_list_syndrome_combined[0]= buffer_pk_in;

		for(int i=0; i<syndrome_kernels; i++){
			pt_list_syndrome_combined[i+1]= buffer_e_in_list[i];
		}





#endif



#ifdef SYND_KERNEL

	//Initalize the buffers/pointers that will be used by all kernels
	buffer_L_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*SYS_N, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_L_in");
		return EXIT_FAILURE;
	}
	#endif

	buffer_r_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_r_in");
		return EXIT_FAILURE;
	}
	#endif



	ptr_L_in = (gf *) clEnqueueMapBuffer(commands, buffer_L_in, true, CL_MAP_WRITE, 0, sizeof(gf)*SYS_N, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_a");
		return EXIT_FAILURE;
	}
	#endif


	ptr_r_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_r_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_r_out");
		return EXIT_FAILURE;
	}
	#endif

	int synd_idx;
	switch (synd_kernels)
	{
    	case 1:
    		synd_idx = 0;
    		break;
	    case 2:
	    	synd_idx = 1;
	    	break;
    	case 4:
    		synd_idx = 3;
    		break;
	    case 8:
	    	synd_idx = 7;
	    	break;
	    case 12:
	    	synd_idx = 15;
	    	break;
	    case 16:
	    	synd_idx = 27;
	    	break;
	    case 18:
	    	synd_idx = 43;
	    	break;
	    default:
	    	break;
	}



	//Iteratively initialize kernels and their private buffers/pointers
	for(int i=0; i<synd_kernels; i++){

		int index = synd_idx +i;

		synd_kernels_list[i] = clCreateKernel(program, synd_kernels_name_list[index], &err);
		#ifdef OCL_API_DEBUG
		if (!synd_kernels_list[i] || err != CL_SUCCESS) {
			printf("Error: Failed to create compute kernel_synd!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		#endif

		buffer_out_out_list[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to create buffer_out_out");
			return EXIT_FAILURE;
		}
		#endif

		buffer_f_in_list[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to create buffer_f_in");
			return EXIT_FAILURE;
		}
		#endif


		err = clSetKernelArg(synd_kernels_list[i], 0, sizeof(cl_mem), &(buffer_out_out_list[i]));
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_out_out");
			return EXIT_FAILURE;
		}
		#endif

		err = clSetKernelArg(synd_kernels_list[i], 1, sizeof(cl_mem), &(buffer_f_in_list[i]));
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_f_in");
			return EXIT_FAILURE;
		}
		#endif

		err = clSetKernelArg(synd_kernels_list[i], 2, sizeof(cl_mem), &buffer_L_in);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_L_in");
			return EXIT_FAILURE;
		}
		#endif


		err = clSetKernelArg(synd_kernels_list[i], 3, sizeof(cl_mem), &buffer_r_in);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_r_in");
			return EXIT_FAILURE;
		}
		#endif


		ptr_out_out_list[i] = (gf *) clEnqueueMapBuffer(commands, buffer_out_out_list[i], true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("ERROR : %d\n", err);
			printf("FAILED to enqueue map buffer_r_out");
			return EXIT_FAILURE;
		}
		#endif

		ptr_f_in_list[i] = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_list[i], true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("ERROR : %d\n", err);
			printf("FAILED to enqueue map buffer_f");
			return EXIT_FAILURE;
		}
		#endif

	}

	//Put buffers/pointers to lists so they can be used by the host function
	pt_list_synd_combined[0]= buffer_L_in;
	pt_list_synd_combined[1]= buffer_r_in;

	for(int i=0; i<synd_kernels; i++){
		pt_list_synd_combined[2+i]= buffer_f_in_list[i];
	}

	for(int i=0; i<synd_kernels; i++){
		pt_list_synd_combined_out[i]= buffer_out_out_list[i];
	}


#endif


    FILE                *fp_req, *fp_rsp;
    int                 ret_val;
    int i;
    unsigned char *ct = 0;
    unsigned char *ss = 0;
    unsigned char *ss1 = 0;
    unsigned char *pk = 0;
    unsigned char *sk = 0;

    struct timeval start_keygen, end_keygen, start_enc, end_enc, start_dec, end_dec;

    for (i=0; i<48; i++)
        entropy_input[i] = i;
    randombytes_init(entropy_input, NULL, 256);

    for (i=0; i<KATNUM; i++)
        randombytes(seed[i], 48);

    fp_req = fopen("kat_kem.req", "w");
    if (!fp_req)
        return KAT_FILE_OPEN_ERROR;

    for (i=0; i<KATNUM; i++) {
        fprintf(fp_req, "count = %d\n", i);
        fprintBstr(fp_req, "seed = ", seed[i], 48);
        fprintf(fp_req, "pk =\n");
        fprintf(fp_req, "sk =\n");
        fprintf(fp_req, "ct =\n");
        fprintf(fp_req, "ss =\n\n");
    }

    fp_rsp = fopen("kat_kem.rsp", "w");
    if (!fp_rsp)
        return KAT_FILE_OPEN_ERROR;

    fprintf(fp_rsp, "# kem/%s\n\n", crypto_kem_PRIMITIVE);

    for (i=0; i<KATNUM; i++) {
        if (!ct) ct = malloc(crypto_kem_CIPHERTEXTBYTES);
        if (!ct) abort();
        if (!ss) ss = malloc(crypto_kem_BYTES);
        if (!ss) abort();
        if (!ss1) ss1 = malloc(crypto_kem_BYTES);
        if (!ss1) abort();
        if (!pk) pk = malloc(crypto_kem_PUBLICKEYBYTES);
        if (!pk) abort();
        if (!sk) sk = malloc(crypto_kem_SECRETKEYBYTES);
        if (!sk) abort();

        randombytes_init(seed[i], NULL, 256);

        fprintf(fp_rsp, "count = %d\n", i);
        fprintBstr(fp_rsp, "seed = ", seed[i], 48);
       
#ifdef TIME_MEASUREMENT
        struct timeval start_keygen, end_keygen;
        gettimeofday(&start_keygen, NULL);
#endif
        ret_val = crypto_kem_keypair(pk, sk);

#ifdef TIME_MEASUREMENT
        gettimeofday(&end_keygen, NULL);
        get_event_time(&start_keygen, &end_keygen, &sum_keygen, &times_keygen);
#endif


        if (ret_val != 0) {
            fprintf(stderr, "crypto_kem_keypair returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "pk = ", pk, crypto_kem_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk = ", sk, crypto_kem_SECRETKEYBYTES);


#ifdef TIME_MEASUREMENT
		struct timeval start_enc, end_enc;
		gettimeofday(&start_enc, NULL);
#endif
		ret_val = crypto_kem_enc(ct, ss, pk);

#ifdef TIME_MEASUREMENT
		gettimeofday(&end_enc, NULL);
		get_event_time(&start_enc, &end_enc, &sum_enc, &times_enc);
#endif


		if (ret_val != 0) {
			fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
			return KAT_CRYPTO_FAILURE;
		}
		fprintBstr(fp_rsp, "ct = ", ct, crypto_kem_CIPHERTEXTBYTES);
		fprintBstr(fp_rsp, "ss = ", ss, crypto_kem_BYTES);

		fprintf(fp_rsp, "\n");

#ifdef TIME_MEASUREMENT
		struct timeval start_dec, end_dec;
		gettimeofday(&start_dec, NULL);
#endif

		ret_val =  crypto_kem_dec(ss1, ct, sk);


#ifdef TIME_MEASUREMENT
		gettimeofday(&end_dec, NULL);
		get_event_time(&start_dec, &end_dec, &sum_dec, &times_dec);
#endif


		if (ret_val != 0) {
			fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
			return KAT_CRYPTO_FAILURE;
		}

		if ( memcmp(ss, ss1, crypto_kem_BYTES) ) {
			fprintf(stderr, "crypto_kem_dec returned bad 'ss' value\n");
			return KAT_CRYPTO_FAILURE;
		}


    }
	

#ifdef TIME_MEASUREMENT
#ifdef GAUSSIAN_ELIMINATION_KERNEL
	printf("\n***************ELIM KERNEL***************\n");
	printf("Kernel execution time\n");
	print_kernel_execution_time(sum_list_elim_kernel, &times_elim, 1);
	printf("To kernel migration time ");
	print_kernel_execution_time(sum_list_elim_tokern, &times_elim_tokern, 1);
	printf("To host migration time ");
	print_kernel_execution_time(sum_list_elim_tohost, &times_elim_tohost, 1);
	printf("Pk loop part ");
	print_event_execution_time(&sum_pk_loop, &times_pk_loop);
	printf("Parallel part ");
	print_event_execution_time(&sum_parallel, &times_parallel);
	printf("While KeyGen part ");
	print_event_execution_time(&sum_while_pk_loop, &times_while_pk_loop);
#endif


#ifdef SYND_KERNEL
	printf("\n***************SYND KERNEL***************\n");
	printf("Kernel execution time\n");
	print_kernel_execution_time(sum_list_synd_kernel, &times_synd, synd_kernels);
	printf("To kernel migration time ");
	print_kernel_execution_time(sum_list_synd_tokern, &times_synd_tokern, 1);
	printf("To host migration time ");
	print_kernel_execution_time(sum_list_synd_tohost, &times_synd_tohost, 1);
	printf("Synd all kernels\n");
	print_event_execution_time(&sum_synd_kernels, &times_synd_kernels);
	printf("Synd host function ");
	print_event_execution_time(&sum_total_synd, &times_total_synd);
#endif

#ifdef SYNDROME_KERNEL
	printf("\n***************SYNDROME KERNEL***************\n");
	printf("Kernel execution time\n");
	print_kernel_execution_time(sum_list_syndrome_kernel, &times_syndrome, syndrome_kernels);

	printf("Kernel execution time\n");
	print_kernel_execution_time(sum_list_syndrome_kernel, &times_syndrome, syndrome_kernels);

	printf("To kernel migration time\n");
	print_kernel_execution_time(sum_list_syndrome_tokern, &times_syndrome_tokern, 1);
	printf("To host migration time\n");
	print_kernel_execution_time(sum_list_syndrome_tohost, &times_syndrome_tohost, 1);
	printf("Syndrome all kernels\n");
	print_event_execution_time(&sum_syndrome_kernels, &times_syndrome_kernels);
	printf("Syndrome host function\n");
	print_event_execution_time(&sum_total_syndrome, &times_total_syndrome);
#endif

#ifdef KEM_PARTS_MEASUREMENT
	printf("\n***************KEM PARTS***************\n");
	printf("Key Generation Part ");
	print_event_execution_time(&sum_keygen, &times_keygen);
	printf("Encapsulate Part ");
	print_event_execution_time(&sum_enc, &times_enc);
	printf("Decapsulate Part ");
	print_event_execution_time(&sum_dec, &times_dec);
	printf("Encrypt Part ");
	print_event_execution_time(&sum_encrypt, &times_encrypt);
	printf("Decrypt Part ");
	print_event_execution_time(&sum_decrypt, &times_decrypt);
#endif

#endif

	#ifdef GAUSSIAN_ELIMINATION_KERNEL
	clReleaseKernel(kernel_gaussian_elimination);
	clEnqueueUnmapMemObject(commands, buffer_mat_in, ptr_mat_in, 0, NULL, NULL);
	clEnqueueUnmapMemObject(commands, buffer_mat_out, NULL, 0, NULL, NULL);
	#endif

	#ifdef SYNDROME_KERNEL
	clReleaseKernel(syndrome_kernels_list[0]);
	clEnqueueUnmapMemObject(commands, buffer_pk_in, ptr_pk_in, 0, NULL, NULL);
	clEnqueueUnmapMemObject(commands, buffer_e_in_list[0], ptr_e_in_list[0], 0, NULL, NULL);
	clEnqueueUnmapMemObject(commands, buffer_s_out, ptr_s_out, 0, NULL, NULL);
	#endif


    clReleaseDevice(device_id);
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return KAT_SUCCESS;
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long i;

	fprintf(fp, "%s", S);

	for ( i=0; i<L; i++ )
		fprintf(fp, "%02X", A[i]);

	if ( L == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}
