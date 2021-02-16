#include "params.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>

gf gf_add_kernel(gf in0, gf in1)
{
	return in0 ^ in1;
}

gf gf_mul_kernel(gf in0, gf in1)
{
	int i;

	uint32_t tmp;
//	, tmp1, tmp2;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;

	t0 = in0;
	t1 = in1;


	tmp = t0 * (t1 & 1);


	for (i = 1; i < GFBITS; i++){
	#pragma HLS unroll factor=4
//#pragma HLS RESOURCE variable=tmp2 core=Mul_lut
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

static inline gf gf_sq_kernel(gf in)
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

gf gf_inv_kernel(gf in)
{
	gf tmp_11;
	gf tmp_1111;

	gf out = in;


	out = gf_sq_kernel(out);
	tmp_11 = gf_mul_kernel(out, in); // 11

	out = gf_sq_kernel(tmp_11);
	out = gf_sq_kernel(out);
	tmp_1111 = gf_mul_kernel(out, tmp_11); // 1111

	out = gf_sq_kernel(tmp_1111);

	out = gf_sq_kernel(out);
	out = gf_sq_kernel(out);
	out = gf_sq_kernel(out);

	out = gf_mul_kernel(out, tmp_1111); // 11111111

	out = gf_sq_kernel(out);
	out = gf_sq_kernel(out);
	out = gf_mul_kernel(out, tmp_11); // 1111111111

	out = gf_sq_kernel(out);
	out = gf_mul_kernel(out, in); // 11111111111

	return gf_sq_kernel(out); // 111111111110
}

gf eval_inner(gf *f, gf a)
{
        int i;
        gf r;

        r = f[ SYS_T ];

        for (i = SYS_T-1; i >= 0; i--)
        {
		#pragma HLS PIPELINE II=3
		#pragma HLS unroll factor=2
                r = gf_mul_kernel(r, a) ^ f[i];
        }

        return r;
}

void synd_kernel(gf *out_out, gf *f_in, gf *L_in, unsigned char *r_in)
{

	#pragma HLS INTERFACE m_axi     port=out_out  offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi     port=f_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=L_in     offset=slave bundle=gmem2
	#pragma HLS INTERFACE m_axi     port=r_in     offset=slave bundle=gmem3
	#pragma HLS INTERFACE s_axilite port=out_out            bundle=control
	#pragma HLS INTERFACE s_axilite port=f_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=L_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=r_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		    bundle=control

	#pragma HLS inline recursive

	int i, j;
	gf e, e_inv, c;
	gf local_out[2*SYS_T];
	gf local_f[SYS_T+1];
	gf local_L[SYS_N/2];
	unsigned char local_r[MAT_COLS/2];

	gf e_mat[SYS_N];

	#pragma HLS ARRAY_PARTITION variable=local_out cyclic factor=2
	#pragma HLS ARRAY_PARTITION variable=local_L cyclic factor=4
	#pragma HLS ARRAY_PARTITION variable=e_mat cyclic factor=4

	//READ into local vars

	LOOP_LOAD_FROM_BRAM_F:for (i=0;i<=SYS_T;i++){
	#pragma HLS PIPELINE II=1
		local_f[i] = *(f_in+i);
	}

	LOOP_LOAD_FROM_BRAM_L:for (i=0;i<SYS_N/2;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=4
		local_L[i] = *(L_in+i);
	}


	LOOP_LOAD_FROM_BRAM_R:for (i=0;i<MAT_COLS/2;i++){
	#pragma HLS PIPELINE II=1
		local_r[i] = *(r_in+i);
	}


	//READ into local vars END
	LOOP_EVAL:
	for(int i=0; i <SYS_N/2; i++){
	#pragma HLS PIPELINE
		e_mat[i] = eval_inner(local_f, local_L[i]);
	}



	LOOP_MAIN_OUTER:
	for (uint i = 0; i < SYS_N/2; i++)
	{
		c = (local_r[i>>3] >> (i%(uint)8)) & 1;
		e_inv = gf_inv_kernel(gf_mul_kernel(e_mat[i],e_mat[i]));

		LOOP_MAIN_INNER:
		for (int j = 0; j < 2*SYS_T; j++)
		{
		#pragma HLS DEPENDENCE inter variable=local_out false
		#pragma HLS PIPELINE II=2
		#pragma HLS unroll factor=2

			if(i==0){
				local_out[j] = gf_mul_kernel(e_inv, c);
			}else{
				local_out[j] ^= gf_mul_kernel(e_inv, c);

			}
			e_inv = gf_mul_kernel(e_inv, local_L[i]);

		}
	}

	//WRITE to local memory

	LOOP_WRITE_TO_BRAM_OUT:for(i=0;i<(2*SYS_T);i++){
	#pragma HLS PIPELINE II=1
		*(out_out+i) = local_out[i];
	}
}
