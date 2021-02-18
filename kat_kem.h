// Global variables to be used by the computational kernels

#include "gf.h"

extern cl_platform_id platform_id;
extern cl_device_id device_id;
extern cl_kernel kernel_gf_add;
extern cl_kernel kernel_gf_mul;
extern cl_kernel kernel_gf_sq;
extern cl_kernel kernel_gaussian_elimination;
extern cl_kernel kernel_eval;
extern cl_kernel kernel_syndrome;
extern cl_kernel kernel_syndrome_2;
extern cl_kernel kernel_syndrome_3;
extern cl_kernel kernel_syndrome_4;
extern cl_kernel kernel_synd;
extern cl_kernel kernel_synd_2;
extern cl_program program;
extern cl_context context;
extern cl_uint DATA_SIZE;
extern cl_int err;
extern cl_command_queue commands;
extern cl_mem buffer_in0;
extern cl_mem buffer_in1;
extern cl_mem buffer_output;


//eval
extern cl_mem buffer_f_in;
extern cl_mem buffer_a_in;
extern cl_mem buffer_r_out;
extern gf *ptr_f_in;
extern gf *ptr_a_in;
extern gf *ptr_r_out;
extern cl_mem pt_list_eval[3];

//syndrome
extern cl_mem buffer_pk_in;
extern cl_mem buffer_e_in;
extern cl_mem buffer_s_out;
extern unsigned char *ptr_pk_in;
extern unsigned char *ptr_e_in;
extern unsigned char *ptr_s_out;
extern cl_mem pt_list_syndrome[3];

extern cl_mem buffer_pk_in_2;
extern cl_mem buffer_e_in_2;
extern cl_mem buffer_s_out_2;
extern unsigned char *ptr_pk_in_2;
extern unsigned char *ptr_e_in_2;
extern unsigned char *ptr_s_out_2;
extern cl_mem pt_list_syndrome_2[3];

extern cl_mem buffer_pk_in_3;
extern cl_mem buffer_e_in_3;
extern cl_mem buffer_s_out_3;
extern unsigned char *ptr_pk_in_3;
extern unsigned char *ptr_e_in_3;
extern unsigned char *ptr_s_out_3;
extern cl_mem pt_list_syndrome_3[3];

extern cl_mem buffer_pk_in_4;
extern cl_mem buffer_e_in_4;
extern cl_mem buffer_s_out_4;
extern unsigned char *ptr_pk_in_4;
extern unsigned char *ptr_e_in_4;
extern unsigned char *ptr_s_out_4;
extern cl_mem pt_list_syndrome_4[3];
extern cl_mem pt_list_syndrome_combined[5];
extern cl_mem pt_list_syndrome_combined_out[4];

//elim
extern cl_mem buffer_mat_in;
extern cl_mem buffer_mat_out;
extern unsigned char *ptr_mat_in;
extern unsigned char *ptr_mat_out;

//synd
extern cl_mem buffer_out_out;
extern cl_mem buffer_f_in;
extern cl_mem buffer_L_in;
extern cl_mem buffer_r_in;
extern gf *ptr_out_out;
extern gf *ptr_f_in;
extern gf *ptr_L_in;
extern unsigned char *ptr_r_in;
extern cl_mem pt_list_synd[4];
extern cl_mem pt_list_synd_combined[5];
extern cl_mem pt_list_synd_combined_out[2];

extern cl_mem buffer_out_out_2;
extern cl_mem buffer_f_in_2;
extern cl_mem buffer_L_in_2;
extern cl_mem buffer_r_in_2;
extern gf *ptr_out_out_2;
extern gf *ptr_f_in_2;
extern gf *ptr_L_in_2;
extern unsigned char *ptr_r_in_2;
extern cl_mem pt_list_synd_2[4];


