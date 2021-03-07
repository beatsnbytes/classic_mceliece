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
#endif

#ifdef EVAL_KERNEL
cl_kernel kernel_eval;
#endif

#ifdef SYNDROME_KERNEL
cl_kernel kernel_syndrome;
cl_kernel kernel_syndrome_2;
cl_kernel kernel_syndrome_3;
cl_kernel kernel_syndrome_4;
#endif

#ifdef SYND_KERNEL
int synd_kernels = 4 ;

cl_kernel kernel_synd_1;
cl_kernel kernel_synd_2;
cl_kernel kernel_synd_3;
cl_kernel kernel_synd_4;
cl_kernel kernel_synd_5;
cl_kernel kernel_synd_6;
cl_kernel kernel_synd_7;
cl_kernel kernel_synd_8;
cl_kernel synd_kernels_list[8];
const char *synd_kernels_name_list[15];

cl_mem buffer_out_out;
cl_mem buffer_out_out_2;
cl_mem buffer_out_out_3;
cl_mem buffer_out_out_4;
cl_mem buffer_out_out_5;
cl_mem buffer_out_out_6;
cl_mem buffer_out_out_7;
cl_mem buffer_out_out_8;
cl_mem buffer_out_out_list[8];

gf *ptr_out_out;
gf *ptr_out_out_2;
gf *ptr_out_out_3;
gf *ptr_out_out_4;
gf *ptr_out_out_5;
gf *ptr_out_out_6;
gf *ptr_out_out_7;
gf *ptr_out_out_8;
gf * ptr_out_out_list[8];

cl_mem buffer_f_in;
cl_mem buffer_f_in_2;
cl_mem buffer_f_in_3;
cl_mem buffer_f_in_4;
cl_mem buffer_f_in_5;
cl_mem buffer_f_in_6;
cl_mem buffer_f_in_7;
cl_mem buffer_f_in_8;
cl_mem buffer_f_in_list[8];

gf *ptr_f_in;
gf *ptr_f_in_2;
gf *ptr_f_in_3;
gf *ptr_f_in_4;
gf *ptr_f_in_5;
gf *ptr_f_in_6;
gf *ptr_f_in_7;
gf *ptr_f_in_8;
gf *ptr_f_in_list[8];

cl_mem buffer_L_in;
gf *ptr_L_in;

cl_mem buffer_r_in;
unsigned char *ptr_r_in;

cl_mem pt_list_synd_combined[10];

cl_mem pt_list_synd_combined_out[8];


#endif

//syndrome
cl_mem buffer_pk_in;
cl_mem buffer_e_in;
cl_mem buffer_s_out;
unsigned char *ptr_pk_in;
unsigned char *ptr_e_in;
unsigned char *ptr_s_out;
cl_mem pt_list_syndrome[3];
cl_mem pt_list_syndrome_combined[2];
cl_mem pt_list_syndrome_combined_out[4];

cl_mem buffer_pk_in_2;
cl_mem buffer_e_in_2;
cl_mem buffer_s_out_2;
unsigned char *ptr_pk_in_2;
unsigned char *ptr_e_in_2;
unsigned char *ptr_s_out_2;
cl_mem pt_list_syndrome_2[3];

cl_mem buffer_pk_in_3;
cl_mem buffer_e_in_3;
cl_mem buffer_s_out_3;
unsigned char *ptr_pk_in_3;
unsigned char *ptr_e_in_3;
unsigned char *ptr_s_out_3;
cl_mem pt_list_syndrome_3[3];

cl_mem buffer_pk_in_4;
cl_mem buffer_e_in_4;
cl_mem buffer_s_out_4;
unsigned char *ptr_pk_in_4;
unsigned char *ptr_e_in_4;
unsigned char *ptr_s_out_4;
cl_mem pt_list_syndrome_4[3];

cl_mem buffer_e_in_5;
cl_mem buffer_e_in_6;
cl_mem buffer_e_in_7;
cl_mem buffer_e_in_8;

unsigned char *ptr_pk_in_5;
unsigned char *ptr_pk_in_6;
unsigned char *ptr_pk_in_7;
unsigned char *ptr_pk_in_8;

cl_mem buffer_pk_in_5;
cl_mem buffer_pk_in_6;
cl_mem buffer_pk_in_7;
cl_mem buffer_pk_in_8;

unsigned char *ptr_e_in_5;
unsigned char *ptr_e_in_6;
unsigned char *ptr_e_in_7;
unsigned char *ptr_e_in_8;


//Gaussian Elimination
cl_mem buffer_mat_in;
cl_mem buffer_mat_out;
unsigned char *ptr_mat_in;
unsigned char *ptr_mat_out;
unsigned int *ptr_fail;
unsigned char *success_info_host_ptr;
cl_buffer_region region_success_info;
cl_mem buffer_success_info;
cl_mem buffer_fail;



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

#ifdef EVAL_KERNEL
	kernel_eval = clCreateKernel(program, "eval_kernel", &err);
	#ifdef OCL_API_DEBUG
	if (!kernel_eval || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel_eval!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif


	buffer_f_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_f_in");
		return EXIT_FAILURE;
	}
	#endif

	buffer_a_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*SYS_N, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_a_in");
		return EXIT_FAILURE;
	}
	#endif

	buffer_r_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*SYS_N, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_r_out");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_eval, 0, sizeof(cl_mem), &buffer_f_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_eval, 1, sizeof(cl_mem), &buffer_a_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_a");
		return EXIT_FAILURE;
	}
	#endif


	err = clSetKernelArg(kernel_eval, 2, sizeof(cl_mem), &buffer_r_out);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_a");
		return EXIT_FAILURE;
	}
	#endif

	ptr_f_in = (gf *) clEnqueueMapBuffer(commands, buffer_f_in, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_f");
		return EXIT_FAILURE;
	}
	#endif

	ptr_a_in = (gf *) clEnqueueMapBuffer(commands, buffer_a_in, true, CL_MAP_WRITE, 0, sizeof(gf)*SYS_N, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_a");
		return EXIT_FAILURE;
	}
	#endif


	ptr_r_out = (gf *) clEnqueueMapBuffer(commands, buffer_r_out, true, CL_MAP_READ, 0, sizeof(gf)*SYS_N, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_r_out");
		return EXIT_FAILURE;
	}
	#endif

	pt_list_eval[0] = buffer_f_in;
	pt_list_eval[1] = buffer_a_in;
	pt_list_eval[2] = buffer_r_out;


#endif


#ifdef SYNDROME_KERNEL
	kernel_syndrome = clCreateKernel(program, "syndrome_kernel", &err);
	#ifdef OCL_API_DEBUG
	if (!kernel_syndrome || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel_syndrome!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif


	buffer_pk_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_pk_in");
		return EXIT_FAILURE;
	}
	#endif

	buffer_e_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_e_in");
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

	err = clSetKernelArg(kernel_syndrome, 0, sizeof(cl_mem), &buffer_pk_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_pk_in");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_syndrome, 1, sizeof(cl_mem), &buffer_e_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_e_in");
		return EXIT_FAILURE;
	}
	#endif


	err = clSetKernelArg(kernel_syndrome, 2, sizeof(cl_mem), &buffer_s_out);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_s_out");
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

	ptr_e_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_e_in");
		return EXIT_FAILURE;
	}
	#endif


	ptr_s_out = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out, true, CL_MAP_READ, 0, sizeof(unsigned char)*SYND_BYTES, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_s_out");
		return EXIT_FAILURE;
	}
	#endif

	pt_list_syndrome[0] = buffer_pk_in;
	pt_list_syndrome[1] = buffer_e_in;
	pt_list_syndrome[2] = buffer_s_out;


#endif


//#2

#ifdef SYNDROME_KERNEL
	kernel_syndrome_2 = clCreateKernel(program, "syndrome_kernel_2", &err);
	#ifdef OCL_API_DEBUG
	if (!kernel_syndrome || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel_syndrome!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif


//	buffer_pk_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_pk_in");
//		return EXIT_FAILURE;
//	}
//	#endif


//	buffer_e_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_e_in");
//		return EXIT_FAILURE;
//	}
//	#endif

//	buffer_s_out_2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*SYND_BYTES/4, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_s_out");
//		return EXIT_FAILURE;
//	}
//	#endif

	err = clSetKernelArg(kernel_syndrome_2, 0, sizeof(cl_mem), &buffer_pk_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_pk_in");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(kernel_syndrome_2, 1, sizeof(cl_mem), &buffer_e_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_e_in");
		return EXIT_FAILURE;
	}
	#endif


	err = clSetKernelArg(kernel_syndrome_2, 2, sizeof(cl_mem), &buffer_s_out);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_s_out");
		return EXIT_FAILURE;
	}
	#endif

//	ptr_pk_in_2 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_pk_in_2, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_pk_in");
//		return EXIT_FAILURE;
//	}
//	#endif

//	ptr_e_in_2 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_in_2, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_e_in");
//		return EXIT_FAILURE;
//	}
//	#endif

//	ptr_s_out_2 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out_2, true, CL_MAP_READ, 0, sizeof(unsigned char)*SYND_BYTES/4, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_s_out");
//		return EXIT_FAILURE;
//	}
//	#endif

	pt_list_syndrome_2[0] = buffer_pk_in;
	pt_list_syndrome_2[1] = buffer_e_in;
	pt_list_syndrome_2[2] = buffer_s_out;


#endif


	//#3

	#ifdef SYNDROME_KERNEL
		kernel_syndrome_3 = clCreateKernel(program, "syndrome_kernel_3", &err);
		#ifdef OCL_API_DEBUG
		if (!kernel_syndrome || err != CL_SUCCESS) {
			printf("Error: Failed to create compute kernel_syndrome!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		#endif


//		buffer_pk_in_3 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("FAILED to create buffer_pk_in");
//			return EXIT_FAILURE;
//		}
//		#endif

//		buffer_e_in_3 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("FAILED to create buffer_e_in");
//			return EXIT_FAILURE;
//		}
//		#endif

//		buffer_s_out_3 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*SYND_BYTES/4, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("FAILED to create buffer_s_out");
//			return EXIT_FAILURE;
//		}
//		#endif

		err = clSetKernelArg(kernel_syndrome_3, 0, sizeof(cl_mem), &buffer_pk_in);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_pk_in");
			return EXIT_FAILURE;
		}
		#endif

		err = clSetKernelArg(kernel_syndrome_3, 1, sizeof(cl_mem), &buffer_e_in);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_e_in");
			return EXIT_FAILURE;
		}
		#endif


		err = clSetKernelArg(kernel_syndrome_3, 2, sizeof(cl_mem), &buffer_s_out);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("FAILED to set kernel arguments for buffer_s_out");
			return EXIT_FAILURE;
		}
		#endif

//		ptr_pk_in_3 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_pk_in_3, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, 0, NULL, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("ERROR : %d\n", err);
//			printf("FAILED to enqueue map buffer_pk_in");
//			return EXIT_FAILURE;
//		}
//		#endif

//		ptr_e_in_3 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_in_3, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("ERROR : %d\n", err);
//			printf("FAILED to enqueue map buffer_e_in");
//			return EXIT_FAILURE;
//		}
//		#endif


//		ptr_s_out_3 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out_3, true, CL_MAP_READ, 0, sizeof(unsigned char)*SYND_BYTES/4, 0, NULL, NULL, &err);
//		#ifdef OCL_API_DEBUG
//		if (err != CL_SUCCESS) {
//			printf("ERROR : %d\n", err);
//			printf("FAILED to enqueue map buffer_s_out");
//			return EXIT_FAILURE;
//		}
//		#endif

		pt_list_syndrome_3[0] = buffer_pk_in;
		pt_list_syndrome_3[1] = buffer_e_in;
		pt_list_syndrome_3[2] = buffer_s_out;


	#endif


		//#4

		#ifdef SYNDROME_KERNEL
			kernel_syndrome_4 = clCreateKernel(program, "syndrome_kernel_4", &err);
			#ifdef OCL_API_DEBUG
			if (!kernel_syndrome || err != CL_SUCCESS) {
				printf("Error: Failed to create compute kernel_syndrome!\n");
				printf("Test failed\n");
				return EXIT_FAILURE;
			}
			#endif


//			buffer_pk_in_4 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("FAILED to create buffer_pk_in");
//				return EXIT_FAILURE;
//			}
//			#endif

//			buffer_e_in_4 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("FAILED to create buffer_e_in");
//				return EXIT_FAILURE;
//			}
//			#endif

//			buffer_s_out_4 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*SYND_BYTES/4, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("FAILED to create buffer_s_out");
//				return EXIT_FAILURE;
//			}
//			#endif

			err = clSetKernelArg(kernel_syndrome_4, 0, sizeof(cl_mem), &buffer_pk_in);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_pk_in");
				return EXIT_FAILURE;
			}
			#endif

			err = clSetKernelArg(kernel_syndrome_4, 1, sizeof(cl_mem), &buffer_e_in);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_e_in");
				return EXIT_FAILURE;
			}
			#endif


			err = clSetKernelArg(kernel_syndrome_4, 2, sizeof(cl_mem), &buffer_s_out);
			#ifdef OCL_API_DEBUG
			if (err != CL_SUCCESS) {
				printf("FAILED to set kernel arguments for buffer_s_out");
				return EXIT_FAILURE;
			}
			#endif

//			ptr_pk_in_4 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_pk_in_4, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES/4, 0, NULL, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("ERROR : %d\n", err);
//				printf("FAILED to enqueue map buffer_pk_in");
//				return EXIT_FAILURE;
//			}
//			#endif

//			ptr_e_in_4 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_in_4, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("ERROR : %d\n", err);
//				printf("FAILED to enqueue map buffer_e_in");
//				return EXIT_FAILURE;
//			}
//			#endif


//			ptr_s_out_4 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out_4, true, CL_MAP_READ, 0, sizeof(unsigned char)*SYND_BYTES/4, 0, NULL, NULL, &err);
//			#ifdef OCL_API_DEBUG
//			if (err != CL_SUCCESS) {
//				printf("ERROR : %d\n", err);
//				printf("FAILED to enqueue map buffer_s_out");
//				return EXIT_FAILURE;
//			}
//			#endif

			pt_list_syndrome_4[0] = buffer_pk_in;
			pt_list_syndrome_4[1] = buffer_e_in;
			pt_list_syndrome_4[2] = buffer_s_out;


			pt_list_syndrome_combined[0]= buffer_pk_in;
//			pt_list_syndrome_combined[1]= buffer_pk_in_2;
//			pt_list_syndrome_combined[2]= buffer_pk_in_3;
//			pt_list_syndrome_combined[3]= buffer_pk_in_4;
			pt_list_syndrome_combined[1]= buffer_e_in;



		#endif




#ifdef SYND_KERNEL

//synd_kernels_list[0] = kernel_synd_1;
//synd_kernels_list[1] = kernel_synd_2;
//synd_kernels_list[2] = kernel_synd_3;
//synd_kernels_list[3] = kernel_synd_4;
//synd_kernels_list[4] = kernel_synd_5;
//synd_kernels_list[5] = kernel_synd_6;
//synd_kernels_list[6] = kernel_synd_7,
//synd_kernels_list[7] = kernel_synd_8;

synd_kernels_name_list[0] = "synd_kernel1_1";
synd_kernels_name_list[1] ="synd_kernel2_1";
synd_kernels_name_list[2]="synd_kernel2_2";
synd_kernels_name_list[3]="synd_kernel4_1";
synd_kernels_name_list[4]="synd_kernel4_2";
synd_kernels_name_list[5]="synd_kernel4_3";
synd_kernels_name_list[6]="synd_kernel4_4";
synd_kernels_name_list[7]="synd_kernel8_1";
synd_kernels_name_list[8]="synd_kernel8_2";
synd_kernels_name_list[9]="synd_kernel8_3";
synd_kernels_name_list[10]="synd_kernel8_4";
synd_kernels_name_list[11]="synd_kernel8_5";
synd_kernels_name_list[12]="synd_kernel8_6";
synd_kernels_name_list[13]="synd_kernel8_7";
synd_kernels_name_list[14]="synd_kernel8_8";


//buffer_out_out_list[0] = buffer_out_out;
//buffer_out_out_list[1] = buffer_out_out_2;
//buffer_out_out_list[2] = buffer_out_out_3;
//buffer_out_out_list[3] = buffer_out_out_4;
//buffer_out_out_list[4] = buffer_out_out_5;
//buffer_out_out_list[5] = buffer_out_out_6;
//buffer_out_out_list[6] = buffer_out_out_7;
//buffer_out_out_list[7] = buffer_out_out_8;


//ptr_out_out_list[0] = ptr_out_out;
//ptr_out_out_list[1] = ptr_out_out_2;
//ptr_out_out_list[2] = ptr_out_out_3;
//ptr_out_out_list[3] = ptr_out_out_4;
//ptr_out_out_list[4] = ptr_out_out_5;
//ptr_out_out_list[5] = ptr_out_out_6;
//ptr_out_out_list[6] = ptr_out_out_7;
//ptr_out_out_list[7] = ptr_out_out_8;



//buffer_f_in_list[0] = buffer_f_in;
//buffer_f_in_list[1] = buffer_f_in_2;
//buffer_f_in_list[2] = buffer_f_in_3;
//buffer_f_in_list[3] = buffer_f_in_4;
//buffer_f_in_list[4] = buffer_f_in_5;
//buffer_f_in_list[5] = buffer_f_in_6;
//buffer_f_in_list[6] = buffer_f_in_7;
//buffer_f_in_list[7] = buffer_f_in_8;
//
//
//ptr_f_in_list[0] = ptr_f_in;
//ptr_f_in_list[1] = ptr_f_in_2;
//ptr_f_in_list[2] = ptr_f_in_3;
//ptr_f_in_list[3] = ptr_f_in_4;
//ptr_f_in_list[4] = ptr_f_in_5;
//ptr_f_in_list[5] = ptr_f_in_6;
//ptr_f_in_list[6] = ptr_f_in_7;
//ptr_f_in_list[7] = ptr_f_in_8;


////if(synd_kernels==1){
////	kernel_synd1_1 = clCreateKernel(program, "synd_kernel1_1", &err);
////	#ifdef OCL_API_DEBUG
////	if (!kernel_synd1_1 || err != CL_SUCCESS) {
////		printf("Error: Failed to create compute kernel_synd!\n");
////		printf("Test failed\n");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	buffer_out_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to create buffer_out_out");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	buffer_f_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to create buffer_f_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	buffer_L_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*SYS_N, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to create buffer_L_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	buffer_r_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to create buffer_r_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	err = clSetKernelArg(kernel_synd1_1, 0, sizeof(cl_mem), &buffer_out_out);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to set kernel arguments for buffer_out_out");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	err = clSetKernelArg(kernel_synd1_1, 1, sizeof(cl_mem), &buffer_f_in);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to set kernel arguments for buffer_f_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	err = clSetKernelArg(kernel_synd1_1, 2, sizeof(cl_mem), &buffer_L_in);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to set kernel arguments for buffer_L_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////
////	err = clSetKernelArg(kernel_synd1_1, 3, sizeof(cl_mem), &buffer_r_in);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to set kernel arguments for buffer_r_in");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	ptr_out_out = (gf *) clEnqueueMapBuffer(commands, buffer_out_out, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("ERROR : %d\n", err);
////		printf("FAILED to enqueue map buffer_r_out");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	ptr_f_in = (gf *) clEnqueueMapBuffer(commands, buffer_f_in, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("ERROR : %d\n", err);
////		printf("FAILED to enqueue map buffer_f");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////	ptr_L_in = (gf *) clEnqueueMapBuffer(commands, buffer_L_in, true, CL_MAP_WRITE, 0, sizeof(gf)*SYS_N, 0, NULL, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("ERROR : %d\n", err);
////		printf("FAILED to enqueue map buffer_a");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////
////	ptr_r_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_r_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("ERROR : %d\n", err);
////		printf("FAILED to enqueue map buffer_r_out");
////		return EXIT_FAILURE;
////	}
////	#endif
////
////
////}else if (synd_kernels==2){
//
//	//2_1
//	kernel_synd2_1 = clCreateKernel(program, "synd_kernel2_1", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd2_1 || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out2_1 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_f_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_L_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*SYS_N, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_r_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd2_1, 0, sizeof(cl_mem), &buffer_out_out2_1);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd2_1, 1, sizeof(cl_mem), &buffer_f_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd2_1, 2, sizeof(cl_mem), &buffer_L_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd2_1, 3, sizeof(cl_mem), &buffer_r_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_out_out2_1 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out2_1, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_f_in = (gf *) clEnqueueMapBuffer(commands, buffer_f_in, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_f");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_L_in = (gf *) clEnqueueMapBuffer(commands, buffer_L_in, true, CL_MAP_WRITE, 0, sizeof(gf)*SYS_N, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_a");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_r_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_r_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//
//	//2_2
//	kernel_synd2_2 = clCreateKernel(program, "synd_kernel2_2", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd2_2 || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out2_2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_f_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd2_2, 0, sizeof(cl_mem), &buffer_out_out2_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd2_2, 1, sizeof(cl_mem), &buffer_f_in_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd2_2, 2, sizeof(cl_mem), &buffer_L_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd2_2, 3, sizeof(cl_mem), &buffer_r_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_out_out2_2 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out2_2, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_f_in_2 = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_2, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_f");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//} else if (synd_kernels==4){

	//4_1
	synd_kernels_list[0] = clCreateKernel(program, "synd_kernel4_1", &err);
	#ifdef OCL_API_DEBUG
	if (!synd_kernels_list[0] || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel_synd!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	buffer_out_out_list[0] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_out_out");
		return EXIT_FAILURE;
	}
	#endif

	buffer_f_in_list[0] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to create buffer_f_in");
		return EXIT_FAILURE;
	}
	#endif

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


	err = clSetKernelArg(synd_kernels_list[0], 0, sizeof(cl_mem), &buffer_out_out_list[0]);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_out_out");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(synd_kernels_list[0], 1, sizeof(cl_mem), &buffer_f_in_list[0]);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_f_in");
		return EXIT_FAILURE;
	}
	#endif

	err = clSetKernelArg(synd_kernels_list[0], 2, sizeof(cl_mem), &buffer_L_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_L_in");
		return EXIT_FAILURE;
	}
	#endif


	err = clSetKernelArg(synd_kernels_list[0], 3, sizeof(cl_mem), &buffer_r_in);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to set kernel arguments for buffer_r_in");
		return EXIT_FAILURE;
	}
	#endif


	ptr_out_out_list[0] = (gf *) clEnqueueMapBuffer(commands, buffer_out_out_list[0], true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_r_out");
		return EXIT_FAILURE;
	}
	#endif

	ptr_f_in_list[0] = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_list[0], true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("ERROR : %d\n", err);
		printf("FAILED to enqueue map buffer_f");
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

	for(int i=1; i<synd_kernels; i++){

		//4_2
		int index = synd_kernels-1+i;
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

//
//	//4_2
//	kernel_synd4_2 = clCreateKernel(program, "synd_kernel4_2", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd4_2 || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out4_2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_f_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_2, 0, sizeof(cl_mem), &buffer_out_out4_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_2, 1, sizeof(cl_mem), &buffer_f_in_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_2, 2, sizeof(cl_mem), &buffer_L_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_2, 3, sizeof(cl_mem), &buffer_r_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_out_out4_2 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out4_2, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_f_in_2 = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_2, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_f");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	//4_3
//	kernel_synd4_3 = clCreateKernel(program, "synd_kernel4_3", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd4_3 || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out4_3 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_f_in_3 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_3, 0, sizeof(cl_mem), &buffer_out_out4_3);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_3, 1, sizeof(cl_mem), &buffer_f_in_3);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_3, 2, sizeof(cl_mem), &buffer_L_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_3, 3, sizeof(cl_mem), &buffer_r_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_out_out4_3 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out4_3, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_f_in_3 = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_3, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_f");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	//4_4
//	kernel_synd4_4 = clCreateKernel(program, "synd_kernel4_4", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd4_4 || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out4_4 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_f_in_4 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_4, 0, sizeof(cl_mem), &buffer_out_out4_4);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_4, 1, sizeof(cl_mem), &buffer_f_in_4);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd4_4, 2, sizeof(cl_mem), &buffer_L_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd4_4, 3, sizeof(cl_mem), &buffer_r_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_out_out4_4 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out4_4, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_f_in_4 = (gf *) clEnqueueMapBuffer(commands, buffer_f_in_4, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_f");
//		return EXIT_FAILURE;
//	}
//	#endif
//

//}
//	pt_list_synd_combined[0]= buffer_f_in;
//	pt_list_synd_combined[1]= buffer_L_in;
//	pt_list_synd_combined[2]= buffer_r_in;
//	pt_list_synd_combined[3]= buffer_f_in_2;
//	pt_list_synd_combined[4]= buffer_f_in_3;
//	pt_list_synd_combined[5]= buffer_f_in_4;
//	pt_list_synd_combined[6]= buffer_f_in_5;
//	pt_list_synd_combined[7]= buffer_f_in_6;
//	pt_list_synd_combined[8]= buffer_f_in_7;
//	pt_list_synd_combined[9]= buffer_f_in_8;
//
//	pt_list_synd_combined_out[0]= buffer_out_out;
//	pt_list_synd_combined_out[1]= buffer_out_out_2;
//	pt_list_synd_combined_out[2]= buffer_out_out_3;
//	pt_list_synd_combined_out[3]= buffer_out_out_4;
//	pt_list_synd_combined_out[4]= buffer_out_out_5;
//	pt_list_synd_combined_out[5]= buffer_out_out_6;
//	pt_list_synd_combined_out[6]= buffer_out_out_7;
//	pt_list_synd_combined_out[7]= buffer_out_out_8;

	pt_list_synd_combined[0]= buffer_f_in_list[0];
	pt_list_synd_combined[1]= buffer_L_in;
	pt_list_synd_combined[2]= buffer_r_in;
	pt_list_synd_combined[3]= buffer_f_in_list[1];
	pt_list_synd_combined[4]= buffer_f_in_list[2];
	pt_list_synd_combined[5]= buffer_f_in_list[3];
	pt_list_synd_combined[6]= buffer_f_in_list[4];
	pt_list_synd_combined[7]= buffer_f_in_list[5];
	pt_list_synd_combined[8]= buffer_f_in_list[6];
	pt_list_synd_combined[9]= buffer_f_in_list[7];

	pt_list_synd_combined_out[0]= buffer_out_out_list[0];
	pt_list_synd_combined_out[1]= buffer_out_out_list[1];
	pt_list_synd_combined_out[2]= buffer_out_out_list[2];
	pt_list_synd_combined_out[3]= buffer_out_out_list[3];
	pt_list_synd_combined_out[4]= buffer_out_out_list[4];
	pt_list_synd_combined_out[5]= buffer_out_out_list[5];
	pt_list_synd_combined_out[6]= buffer_out_out_list[6];
	pt_list_synd_combined_out[7]= buffer_out_out_list[7];



#endif


////#2
//
//#ifdef SYND_KERNEL
//	kernel_synd_2 = clCreateKernel(program, "synd_kernel_2", &err);
//	#ifdef OCL_API_DEBUG
//	if (!kernel_synd || err != CL_SUCCESS) {
//		printf("Error: Failed to create compute kernel_synd!\n");
//		printf("Test failed\n");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_out_out_2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(gf)*2*SYS_T, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
////	buffer_f_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*(SYS_T+1), NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("FAILED to create buffer_f_in");
////		return EXIT_FAILURE;
////	}
////	#endif
//
//	buffer_L_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(gf)*SYS_N/2, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	buffer_r_in_2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*MAT_COLS/2, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to create buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd_2, 0, sizeof(cl_mem), &buffer_out_out_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_out_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd_2, 1, sizeof(cl_mem), &buffer_f_in);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_f_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	err = clSetKernelArg(kernel_synd_2, 2, sizeof(cl_mem), &buffer_L_in_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_L_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	err = clSetKernelArg(kernel_synd_2, 3, sizeof(cl_mem), &buffer_r_in_2);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("FAILED to set kernel arguments for buffer_r_in");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	ptr_out_out_2 = (gf *) clEnqueueMapBuffer(commands, buffer_out_out_2, true, CL_MAP_READ, 0, sizeof(gf)*2*SYS_T, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
////	ptr_f_in = (gf *) clEnqueueMapBuffer(commands, buffer_f_in, true, CL_MAP_WRITE, 0, sizeof(gf)*(SYS_T+1), 0, NULL, NULL, &err);
////	#ifdef OCL_API_DEBUG
////	if (err != CL_SUCCESS) {
////		printf("ERROR : %d\n", err);
////		printf("FAILED to enqueue map buffer_f");
////		return EXIT_FAILURE;
////	}
////	#endif
//
//	ptr_L_in_2 = (gf *) clEnqueueMapBuffer(commands, buffer_L_in_2, true, CL_MAP_WRITE, 0, sizeof(gf)*SYS_N/2, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_a");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//
//	ptr_r_in_2 = (unsigned char *) clEnqueueMapBuffer(commands, buffer_r_in_2, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*MAT_COLS/2, 0, NULL, NULL, &err);
//	#ifdef OCL_API_DEBUG
//	if (err != CL_SUCCESS) {
//		printf("ERROR : %d\n", err);
//		printf("FAILED to enqueue map buffer_r_out");
//		return EXIT_FAILURE;
//	}
//	#endif
//
//	pt_list_synd_2[0] = buffer_out_out_2;
//	pt_list_synd_2[1] = buffer_f_in;
//	pt_list_synd_2[2] = buffer_L_in_2;
//	pt_list_synd_2[3] = buffer_r_in_2;
//
////	cl_mem pt_list_synd_combined[5];
//	pt_list_synd_combined[0]=buffer_f_in;
//	pt_list_synd_combined[1]=buffer_L_in;
//	pt_list_synd_combined[2]=buffer_r_in;
//	pt_list_synd_combined[3]=buffer_L_in_2;
//	pt_list_synd_combined[4]=buffer_r_in_2;
//
////	cl_mem pt_list_synd_combined_out[2];
//	pt_list_synd_combined_out[0] = buffer_out_out;
//	pt_list_synd_combined_out[1] = buffer_out_out_2;
//
//
//#endif






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
        
        //
//         for(int t=0; t<100; t++){
        //
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

//        }//t

    }
	

#ifdef TIME_MEASUREMENT
//	printf("\n***************TIMING RESULTS***************\n");
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
	printf("XOR reduction\n");
	print_event_execution_time(&sum_xor, &times_xor);
	printf("Synd all kernels\n");
	print_event_execution_time(&sum_synd_kernels, &times_synd_kernels);
	printf("Synd host function ");
	print_event_execution_time(&sum_total_synd, &times_total_synd);
#endif

	#ifdef SYNDROME_KERNEL
    printf("Syndrome kernel :Avg Execution time is: %0.3f miliseconds \n", (sum_syndrome/1000000.0)/times_syndrome);
    printf("Syndrome_2 kernel :Avg Execution time is: %0.3f miliseconds \n", (sum_syndrome_2/1000000.0)/times_syndrome_2);
    printf("Syndrome_3 kernel :Avg Execution time is: %0.3f miliseconds \n", (sum_syndrome_3/1000000.0)/times_syndrome_3);
    printf("Syndrome_4 kernel :Avg Execution time is: %0.3f miliseconds \n", (sum_syndrome_4/1000000.0)/times_syndrome_4);
    printf("All Syndrome kernels :Avg Execution time is: %0.3f miliseconds \n", ((sum_syndrome + sum_syndrome_2 + sum_syndrome_3 + sum_syndrome_4)/4000000.0)/(times_syndrome_4));
//    printf("All Syndrome kernels :Avg Execution time is: %0.3f miliseconds \n", ((sum_syndrome + sum_syndrome_2)/2000000.0)/(times_syndrome_2));
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
	clReleaseKernel(kernel_syndrome);
	clEnqueueUnmapMemObject(commands, buffer_pk_in, ptr_pk_in, 0, NULL, NULL);
	clEnqueueUnmapMemObject(commands, buffer_e_in, ptr_e_in, 0, NULL, NULL);
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
