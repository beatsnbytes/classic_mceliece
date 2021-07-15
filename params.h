#ifndef PARAMS_H
#define PARAMS_H

#define GFBITS 12
#define SYS_N 3488
#define SYS_T 64
#define MAT_SIZE GFBITS*SYS_T*(SYS_N/8)
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

//#define DATAFLOW_OPT
#undef DATAFLOW_OPT

#define PACK_FACTOR_PK 16
#define PACK_BITWIDTH_PK (PACK_FACTOR_PK*8)
#define PACK_FACTOR_E 4
#define PACK_BITWIDTH_E (PACK_FACTOR_E*8)

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

//#define SYNDROME_KERNEL
#undef SYNDROME_KERNEL

//#define SYND_KERNEL
#undef SYND_KERNEL

#define ROOT_KERNEL
//#undef ROOT_KERNEL


typedef struct {
   unsigned char packed_values[PACK_FACTOR_PK];
   }data_packed_pk;

typedef struct {
  unsigned char packed_values[PACK_FACTOR_E];
  }data_packed_e;


#endif

