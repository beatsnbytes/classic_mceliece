#include "params.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>


gf gf_add_kernel_root(gf in0, gf in1)
{
	return in0 ^ in1;
}

gf gf_mul_kernel_root(gf in0, gf in1)
{
	int i;

	uint32_t tmp;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;

	t0 = in0;
	t1 = in1;


	tmp = t0 * (t1 & 1);

	LOOP_MULT:
	for (uint i = 1; i < GFBITS; i++){
	#pragma HLS pipeline

		tmp ^= (t0 * (t1 & (1 << i)));
	}


	t = tmp & 0x7FC000;
	tmp ^= t >> 9;
	tmp ^= t >> 12;

	t = tmp & 0x3000;
	tmp ^= t >> 9;
	tmp ^= t >> 12;

	return tmp & ((1 << GFBITS)-1);
}



void root_kernel1_1(gf* out_out, gf* fr_in, gf* Lr_in)
{

	#pragma HLS INTERFACE m_axi     port=out_out   offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=fr_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=Lr_in     offset=slave bundle=gmem2
	#pragma HLS INTERFACE s_axilite port=out_out             bundle=control
	#pragma HLS INTERFACE s_axilite port=fr_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=Lr_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

	#pragma HLS inline recursive


	int i, j;
	gf r;
	gf local_fr[SYS_T+1];
	gf local_Lr[SYS_N];
	gf local_out[SYS_N];

	#pragma HLS ARRAY_PARTITION variable=local_fr cyclic factor=65
	#pragma HLS ARRAY_PARTITION variable=local_Lr cyclic factor=109
	#pragma HLS ARRAY_PARTITION variable=local_out cyclic factor=109

	LOOP_LOAD_FROM_BRAM_FR:for (int i=0;i<=SYS_T;i++){
	#pragma HLS PIPELINE II=1
//	#pragma HLS unroll factor=2
		local_fr[i] = *(fr_in+i);
	}

	LOOP_LOAD_FROM_BRAM_LR:for (int i=0;i<SYS_N;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_Lr[i] = *(Lr_in+i);
	}

	LOOP_MAIN_OUTTER:
	for (j = 0; j < SYS_N; j++){
	#pragma HLS PIPELINE II=1

		r = local_fr[ SYS_T ];

		LOOP_MAIN_INNER:
		for (i = SYS_T-1; i >= 0; i--)
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=64
		{
			r = gf_mul_kernel_root(r, local_Lr[j]);
			r = gf_add_kernel_root(r, local_fr[i]);
		}

		local_out[j] = r;
	}

	LOOP_WRITE_TO_BRAM_OUT:for(int i=0;i<(SYS_N);i++){
	#pragma HLS PIPELINE II=1
		*(out_out+i) = local_out[i];
	}


}
