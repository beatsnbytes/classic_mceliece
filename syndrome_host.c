

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#include "syndrome_host.h"
#include "params.h"
#include "kat_kem.h"
#include "encrypt.h"
#include "custom_util.h"
#include "sys/time.h"
#include "crypto_kem.h"
#include "params.h"

#ifdef SYNDROME_KERNEL
int syndrome_host(unsigned char *s, const unsigned char *pk, unsigned char  *e)
{

	cl_event events_enq[8], event_migr_tokern, events_migr_tohost;


//	ap_uint<PACK_BITWIDTH_PK> zeros = 0;
//	ap_uint<PACK_BITWIDTH_PK> buff_array[PACK_FACTOR_PK];
//	for(int i=0; i<MAT_ROWS; i++){
//	   for(int j=0; j<(PK_ROW_BYTES/PACK_FACTOR_PK); j++){
//
//		   *(ptr_pk_in+i*(PK_ROW_BYTES/PACK_FACTOR_PK)+j) = zeros;
//		   for(int k=0; k<PACK_FACTOR_PK; k++){
//			   buff_array[k] = zeros;
//			   buff_array[k] = *(pk+i*PK_ROW_BYTES+(j*PACK_FACTOR_PK)+k);
//			   *(ptr_pk_in+i*(PK_ROW_BYTES/PACK_FACTOR_PK)+j) |= (buff_array[k] << (8*k));
//		   }
//	   }
//	}


//for(int i=0; i<MAT_ROWS; i++){
//   for(int j=0; j<(PK_ROW_BYTES/PACK_FACTOR_PK); j++){
//	   for(int k=0; k<PACK_FACTOR_PK; k++){
//	   (*(ptr_pk_in + i*(PK_ROW_BYTES/PACK_FACTOR_PK)+ j)).packed_pk_values[k] = *(pk+i*PK_ROW_BYTES+(j*PACK_FACTOR_PK)+k);
//	   }
//   }
//}

//seems that by puting discrete variables and not matrices works
memcpy(ptr_pk_in, pk, sizeof(unsigned char)*crypto_kem_PUBLICKEYBYTES);





	for(int i=0; i<syndrome_kernels; i++){
		memcpy(ptr_e_in_list[i], e, sizeof(unsigned char)*(MAT_COLS));
	}





	err = clEnqueueMigrateMemObjects(commands, syndrome_kernels+1, &pt_list_syndrome_combined[0], 0, 0, NULL, &event_migr_tokern);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue pt_list_syndrome\n");
		return EXIT_FAILURE;
	}
	#endif

//	#ifdef TIME_MEASUREMENT
//		clWaitForEvents(1, &event_migr_tokern);
//		struct timeval start_kernel, end_kernel;
//		gettimeofday(&start_kernel, NULL);
//	#endif

	for (int i=0; i<syndrome_kernels; i++){
		err = clEnqueueTask(commands, syndrome_kernels_list[i], 1, &event_migr_tokern, &events_enq[i]);
	}

//	#ifdef TIME_MEASUREMENT
//		clWaitForEvents(syndrome_kernels, &events_enq[0]);
//		gettimeofday(&end_kernel, NULL);
//		get_event_time(&start_kernel, &end_kernel, &sum_syndrome_kernels, &times_syndrome_kernels);
//	#endif

	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to execute kernel\n");
		return EXIT_FAILURE;
	}
	#endif


	err = clEnqueueMigrateMemObjects(commands, (cl_uint)1, &buffer_s_out, CL_MIGRATE_MEM_OBJECT_HOST, syndrome_kernels, events_enq, &events_migr_tohost); //without[0]
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("FAILED to enqueue bufer_res\n");
		return EXIT_FAILURE;
	}
	#endif


	clWaitForEvents(1, &events_migr_tohost);

	memcpy(s, ptr_s_out, sizeof(unsigned char)*SYND_BYTES);



//	unsigned char validate_mat[SYND_BYTES];
//	syndrome_sw_host(validate_mat, pk, e);
//	for (int i=0;i<96;i++){
//		if (validate_mat[i] != *(s+i)){
//			printf("\nERROR in %d: Expected %d, got %d\n", i, validate_mat[i], *(s+i));
//		}
//	}
//	#endif



//#ifdef TIME_MEASUREMENT
//	cl_profile_print(&event_migr_tokern, 1, sum_list_syndrome_tokern, &times_syndrome_tokern);
//	cl_profile_print(&events_enq[0], syndrome_kernels, sum_list_syndrome_kernel, &times_syndrome);
//	cl_profile_print(&events_migr_tohost, 1, sum_list_syndrome_tohost, &times_syndrome_tohost);
//#endif

return 0;
}
#endif


static void syndrome_sw_host(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	unsigned char b, row[SYS_N/8];
	const unsigned char *pk_ptr = pk;

	int i, j;

	for (i = 0; i < SYND_BYTES; i++)
		s[i] = 0;

	for (i = 0; i < PK_NROWS; i++)
	{
		for (j = 0; j < SYS_N/8; j++)
			row[j] = 0;

		for (j = 0; j < PK_ROW_BYTES; j++)
			row[ SYS_N/8 - PK_ROW_BYTES + j ] = pk_ptr[j];

		row[i/8] |= 1 << (i%8);

		b = 0;
		for (j = 0; j < SYS_N/8; j++)
			b ^= row[j] & e[j];

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

		s[ i/8 ] |= (b << (i%8));

		pk_ptr += PK_ROW_BYTES;
	}
}

