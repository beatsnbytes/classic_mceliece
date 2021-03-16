#ifndef PARAMS_H
#define PARAMS_H

#define GFBITS 13
#define SYS_N 6960
#define SYS_T 119

#define MAT_ROWS (SYS_T*GFBITS)
#define MAT_COLS (SYS_N/8)
#define COND_BYTES ((1 << (GFBITS-4))*(2*GFBITS - 1))
#define IRR_BYTES (SYS_T * 2)

#define PK_NROWS (SYS_T*GFBITS) 
#define PK_NCOLS (SYS_N - PK_NROWS)
#define PK_ROW_BYTES ((PK_NCOLS + 7)/8)

#define SYND_BYTES ((PK_NROWS + 7)/8)

#define GFMASK ((1 << GFBITS) - 1)

//Custom definitions

//#define FUNC_CORRECTNESS
#undef FUNC_CORRECTNESS

#define TIME_MEASUREMENT
//#undef TIME_MEASUREMENT

//#define ITER_PRINT
#undef ITERATION_PRINT

#define KEM_PARTS_MEASUREMENT
//#undef KEM_PARTS_MEASUREMENT

#define OCL_API_DEBUG
//#undef OCL_API_DEBUG

//#define GAUSSIAN_ELIMINATION_KERNEL
#undef GAUSSIAN_ELIMINATION_KERNEL

#define SYNDROME_KERNEL
//#undef SYNDROME_KERNEL

//#define SYND_KERNEL
#undef SYND_KERNEL


#endif

