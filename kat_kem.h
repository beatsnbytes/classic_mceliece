// Global variables to be used by the computational kernels

#include "gf.h"
//#include "ap_int.h"

extern cl_platform_id platform_id;
extern cl_device_id device_id;
extern cl_kernel kernel_gaussian_elimination;
extern cl_kernel kernel_eval;
extern cl_kernel kernel_syndrome;
extern cl_kernel kernel_syndrome_2;
extern cl_kernel kernel_syndrome_3;
extern cl_kernel kernel_syndrome_4;
extern cl_kernel kernel_syndrome_5;
extern cl_kernel kernel_syndrome_6;
extern cl_kernel kernel_syndrome_7;
extern cl_kernel kernel_syndrome_8;


extern cl_program program;
extern cl_context context;
extern cl_uint DATA_SIZE;
extern cl_int err;
extern cl_command_queue commands;

//syndrome


extern int syndrome_kernels;
extern cl_mem buffer_pk_in;
extern cl_mem buffer_e_in;
extern cl_mem buffer_s_out;


//extern data_packed_pk *ptr_pk_in;
extern unsigned char *ptr_pk_in;


extern data_packed_e *ptr_e_in;
extern unsigned char *ptr_s_out;
extern cl_mem pt_list_syndrome[3];

//
extern unsigned char *ptr_pk_out;
extern cl_mem buffer_pk_out;
//



extern cl_kernel syndrome_kernels_list[8];

extern cl_mem pt_list_syndrome_combined[9];

extern cl_mem buffer_e_in_list[8];
//extern data_packed_e *ptr_e_in_list[8];
extern unsigned char *ptr_e_in_list[8];

extern cl_mem pt_list_syndrome_combined[9];
extern cl_mem pt_list_syndrome_combined_out[8];

//elim
extern cl_mem buffer_mat_in;
extern cl_mem buffer_mat_out;
extern cl_mem buffer_success_info;
extern cl_mem buffer_fail;
extern unsigned char *ptr_mat_in;
extern unsigned char *ptr_mat_out;
extern unsigned char *success_info_host_ptr;
extern unsigned int *ptr_fail;

//SYND
extern int synd_kernels;
extern cl_kernel synd_kernels_list[11];
extern cl_mem buffer_out_out_list[11];
extern gf * ptr_out_out_list[11];
extern cl_mem buffer_f_in_list[11];
extern gf *ptr_f_in_list[11];
extern cl_mem buffer_L_in;
extern gf *ptr_L_in;
extern cl_mem buffer_r_in;
extern unsigned char *ptr_r_in;
extern cl_mem pt_list_synd_combined[13];
extern cl_mem pt_list_synd_combined_out[11];

//ROOT
extern int root_kernels;
extern cl_kernel root_kernels_list[1];

extern cl_mem buffer_out_out;
extern cl_mem buffer_fr_in;
extern cl_mem buffer_Lr_in;

extern gf * ptr_out_out;
extern gf * ptr_fr_in;
extern gf * ptr_Lr_in;

//BM
extern int bm_kernels;
extern cl_kernel bm_kernels_list[1];

extern cl_mem buffer_outbm_out;
extern cl_mem buffer_sbm_in;

extern gf * ptr_outbm_out;
extern gf * ptr_sbm_in;

//SUPPORT

extern int support_kernels;
extern cl_kernel support_kernels_list[1];

extern cl_mem buffer_ssupp_out;
extern cl_mem buffer_csupp_in;


extern gf * ptr_ssupp_out;
extern unsigned char * ptr_csupp_in;


//INIT
extern int init_kernels;

extern cl_kernel init_kernels_list[1];

extern cl_mem buffer_ginit_out;
extern cl_mem buffer_rinit_out;
extern cl_mem buffer_cinit_in;
extern cl_mem buffer_skinit_in;


extern gf * ptr_ginit_out;
extern unsigned char * ptr_rinit_out;
extern unsigned char * ptr_cinit_in;
extern unsigned char * ptr_skinit_in;

//MID

extern int mid_kernels;
extern cl_kernel mid_kernels_list[1];

extern cl_mem buffer_imagesmid_in;
extern cl_mem buffer_wmid_out;
extern cl_mem buffer_emid_out;


extern gf * ptr_imagesmid_in;
extern int * ptr_wmid_out;
extern unsigned char * ptr_emid_out;

//CHECK

extern int check_kernels;

extern cl_kernel check_kernels_list[1];

extern cl_mem buffer_check_in;
extern cl_mem buffer_wcheck_in;
extern cl_mem buffer_scheck_in;
extern cl_mem buffer_s_cmpcheck_in;


extern uint16_t * ptr_check_out;
extern int * ptr_wcheck_in;
extern gf * ptr_scheck_in;
extern gf * ptr_s_cmpcheck_in;


