#include "params.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>
//#include "ap_cint.h"

gf gf_add_kernel8_5(gf in0, gf in1)
{
	return in0 ^ in1;
}

gf gf_mul_kernel8_5(gf in0, gf in1)
{
	int i;

	uint32_t tmp;
	uint32_t tmp_bit, tmp_mul;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;

	t0 = in0;
	t1 = in1;


	tmp = t0 * (t1 & 1);

	for (uint i = 1; i < GFBITS; i++){
//	#pragma HLS DEPENDENCE inter variable=tmp_mul RAW false
//	#pragma HLS pipeline
//	#pragma HLS unroll
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

static inline gf gf_sq_kernel8_5(gf in)
{
	const uint32_t B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};

	uint32_t x = in;
	uint32_t t;

	x = (x | (x << 8)) & B[3];
	x = (x | (x << 4)) & B[2];
	x = (x | (x << 2)) & B[1];
	x = (x | (x << 1)) & B[0];

	t = x & 0x7FC000;
	x ^= t >> 9;
	x ^= t >> 12;

	t = x & 0x3000;
	x ^= t >> 9;
	x ^= t >> 12;

	return x & ((1 << GFBITS)-1);
}

gf gf_inv_kernel8_5(gf in)
{
	gf tmp_11;
	gf tmp_1111;

	gf out = in;


	out = gf_sq_kernel8_5(out);
	tmp_11 = gf_mul_kernel8_5(out, in); // 11

	out = gf_sq_kernel8_5(tmp_11);
	out = gf_sq_kernel8_5(out);
	tmp_1111 = gf_mul_kernel8_5(out, tmp_11); // 1111

	out = gf_sq_kernel8_5(tmp_1111);

	out = gf_sq_kernel8_5(out);
	out = gf_sq_kernel8_5(out);
	out = gf_sq_kernel8_5(out);

	out = gf_mul_kernel8_5(out, tmp_1111); // 11111111

	out = gf_sq_kernel8_5(out);
	out = gf_sq_kernel8_5(out);
	out = gf_mul_kernel8_5(out, tmp_11); // 1111111111

	out = gf_sq_kernel8_5(out);
	out = gf_mul_kernel8_5(out, in); // 11111111111

	return gf_sq_kernel8_5(out); // 111111111110
}

gf eval_inner8_5(gf *f, gf a)
{
        int i;
        gf r;

        r = f[ SYS_T ];

        for (i = SYS_T-1; i >= 0; i--)
        {
//		#pragma HLS PIPELINE II=1
//		#pragma HLS unroll factor=2
                r = gf_mul_kernel8_5(r, a) ^ f[i];
//                r = gf_add(r, f[i]);
        }

        return r;
}

void synd_kernel8_5(gf *out_out, gf *f_in, gf *L_in, unsigned char *r_in)
{

	#pragma HLS INTERFACE m_axi     port=out_out  offset=slave bundle=gmem5
	#pragma HLS INTERFACE m_axi     port=f_in     offset=slave bundle=gmem5
	#pragma HLS INTERFACE m_axi     port=L_in     offset=slave bundle=gmem5
	#pragma HLS INTERFACE m_axi     port=r_in     offset=slave bundle=gmem5
	#pragma HLS INTERFACE s_axilite port=out_out            bundle=control
	#pragma HLS INTERFACE s_axilite port=f_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=L_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=r_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		    bundle=control

	#pragma HLS inline recursive

	int i, j;
	gf e, e_inv, c;
	gf local_out[5*SYS_T];
	gf local_f[SYS_T+1];
	gf local_L[SYS_N];
	gf tmp_mul_1, tmp_mul_2;
	unsigned char local_r[MAT_COLS];

	gf e_mat[SYS_N];

	#pragma HLS ARRAY_PARTITION variable=local_out cyclic factor=4 //64
	#pragma HLS ARRAY_PARTITION variable=local_L cyclic factor=4 //8
	#pragma HLS ARRAY_PARTITION variable=e_mat cyclic factor=4 //8

	//READ into local vars

	LOOP_LOAD_FROM_BRAM_F:for (uint i=0;i<=SYS_T;i++){
	#pragma HLS PIPELINE II=1
		local_f[i] = *(f_in+i);
	}

	LOOP_LOAD_FROM_BRAM_L:for (uint i=4*SYS_N/8;i<5*SYS_N/8;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=4
		local_L[i] = *(L_in+i);
	}

	LOOP_LOAD_FROM_BRAM_R:for (uint i=4*MAT_COLS/8;i<5*MAT_COLS/8;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_r[i] = *(r_in+i);
	}

	//READ into local vars END
	LOOP_EVAL:
	for(uint i=4*SYS_N/8; i <5*SYS_N/8; i++){
//	#pragma HLS PIPELINE
//	#pragma HLS unroll factor=2
		e_mat[i] = eval_inner8_5(local_f, local_L[i]);
	}



	LOOP_MAIN_OUTER:
	for (uint i = 4*SYS_N/8; i < 5*SYS_N/8; i++) //12
	{
	// #pragma HLS pipeline

		c = (local_r[i>>3] >> (i%8)) & 1;
		e_inv = gf_inv_kernel8_5(gf_mul_kernel8_5(e_mat[i],e_mat[i]));

		LOOP_MAIN_INNER:
		for (uint j = 0; j < 5*SYS_T; j++)//8
		{
//		#pragma HLS DEPENDENCE inter variable=local_out false
		#pragma HLS PIPELINE
//		#pragma HLS unroll factor=32

			if(i==4*SYS_N/8){
				local_out[j] = gf_mul_kernel8_5(e_inv, c);
			}else{
				local_out[j] ^= gf_mul_kernel8_5(e_inv, c);

			}
			e_inv = gf_mul_kernel8_5(e_inv, local_L[i]);

		}
	}

	//WRITE to local memory

	LOOP_WRITE_TO_BRAM_OUT:for(uint i=0;i<(5*SYS_T);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		*(out_out+i) = local_out[i];
	}
}
