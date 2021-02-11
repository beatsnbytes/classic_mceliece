#include "root.h"
#include "params.h"
#include <stdlib.h>
#include <string.h>


void eval_kernel(gf *f_in, gf *a_in, gf *r_out)
{
#pragma HLS INTERFACE m_axi     port=f_in   offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi     port=a_in   offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi     port=r_out  offset=slave bundle=gmem2
#pragma HLS INTERFACE s_axilite port=f_in               bundle=control
#pragma HLS INTERFACE s_axilite port=a_in               bundle=control
#pragma HLS INTERFACE s_axilite port=r_out              bundle=control
#pragma HLS INTERFACE s_axilite port=return 		    bundle=control


	int i, j, k;
	gf r_tmp, r_tmp1;
	gf local_f[SYS_T+1];
	gf local_a[SYS_N];
	gf local_r[SYS_N];
	uint32_t tmp;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;



//#pragma HLS ARRAY_PARTITION variable=local_f cyclic factor=3
#pragma HLS ARRAY_PARTITION variable=local_a cyclic
//#pragma HLS ARRAY_PARTITION variable=local_r cyclic factor=3

	LOOP_LOAD_FROM_BRAM_F:for (i=0;i<=SYS_T;i++){
	#pragma HLS PIPELINE II=1
		local_f[i] = *(f_in+i);
	}

	LOOP_LOAD_FROM_BRAM_A:for (i=0;i<SYS_N;i++){
	#pragma HLS PIPELINE II=1
//		#pragma HLS unroll factor=2
		local_a[i] = *(a_in+i);
	}


	SYS_N_ITERS_LOOP:for (k = 0; k < SYS_N; k++)
	{
		#pragma HLS PIPELINE
		#pragma HLS unroll factor=3

		r_tmp1 = local_f[SYS_T];
		OUTER_LOOP:for (i = SYS_T-1; i >= 0; i--)
		{
			#pragma HLS dependence variable=local_r inter false
//			#pragma HLS unroll factor=2
//			#pragma HLS PIPELINE II=2

			t0 = r_tmp1;
			t1 = local_a[k];

			tmp = t0 * (t1 & 1);

			INNER_MUL_LOOP:for (j = 1; j < GFBITS; j++)
			{
				#pragma HLS unroll
				#pragma HLS PIPELINE II=1
				tmp ^= (t0 * (t1 & (1 << j)));
			}

			t = tmp & 0x7FC000;
			tmp ^= t >> 9;
			tmp ^= t >> 12;

			t = tmp & 0x3000;
			tmp ^= t >> 9;
			tmp ^= t >> 12;

			r_tmp = tmp & ((1 << GFBITS)-1);

			r_tmp1 = r_tmp ^ local_f[i];
		}
		local_r[k] = r_tmp1;
	}

	LOOP_WRITE_TO_BRAM_R:for (i=0;i<SYS_N;i++){
		#pragma HLS unroll factor=2
		*(r_out+i) = local_r[i];
	}

}
