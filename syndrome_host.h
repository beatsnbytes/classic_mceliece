
#ifdef __cplusplus
extern "C" {
#endif

extern int syndrome_host(unsigned char *, const unsigned char *, unsigned char  *);

extern double sum_list_syndrome_tokern[1];
extern double sum_list_syndrome_tohost[1];
extern double sum_list_syndrome_kernel[8];
extern int times_syndrome;
extern int times_syndrome_tohost;
extern int times_syndrome_tokern;

extern double sum_syndrome_kernels;
extern int times_syndrome_kernels;

extern double sum_total_syndrome;
extern int times_total_syndrome;

#ifdef __cplusplus
}
#endif
