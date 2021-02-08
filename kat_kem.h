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
extern cl_kernel kernel_composeinv;
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

extern cl_mem buffer_pk_in;
extern cl_mem buffer_e_in;
extern cl_mem buffer_s_out;
extern unsigned char *ptr_pk_in;
extern unsigned char *ptr_e_in;
extern unsigned char *ptr_s_out;
extern cl_mem pt_list_syndrome[3];

//elim
extern cl_mem buffer_mat_in;
extern cl_mem buffer_mat_out;
extern unsigned char *ptr_mat_in;
extern unsigned char *ptr_mat_out;

//composeinv
extern cl_mem buffer_y_in;
extern cl_mem buffer_x_in;
extern cl_mem buffer_pi_in;
extern cl_mem buffer_mats_out;
extern uint32_t *ptr_y_in;
extern uint32_t *ptr_x_in;
extern uint32_t *ptr_pi_in;
extern uint32_t *ptr_mats_out;
extern cl_mem pt_list_composeinv[4];

