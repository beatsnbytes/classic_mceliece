#include "params.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>
//#include "ap_cint.h"

gf gf_add_kernel29_2(gf in0, gf in1)
{
	return in0 ^ in1;
}

gf gf_mul_kernel29_2(gf in0, gf in1)
{
	int i;

	uint64_t tmp;
	uint64_t t0;
	uint64_t t1;
	uint64_t t;

	t0 = in0;
	t1 = in1;

	tmp = t0 * (t1 & 1);

	for (uint i = 1; i < GFBITS; i++){
//#pragma HLS dependence variable=tmp_mul_mat RAW true
//	#pragma HLS pipeline
//	#pragma HLS unroll
		tmp ^= (t0 * (t1 & (1 << i)));
	}


	t = tmp & 0x1FF0000;
	tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);

	t = tmp & 0x000E000;
	tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);

	return tmp & GFMASK;
}


/* input: field element in */
/* return: (in^2)^2 */
static inline gf gf_sq2_kernel29_2(gf in)
{
	int i;

	const uint64_t B[] = {0x1111111111111111,
	                      0x0303030303030303,
	                      0x000F000F000F000F,
	                      0x000000FF000000FF};

	const uint64_t M[] = {0x0001FF0000000000,
	                      0x000000FF80000000,
	                      0x000000007FC00000,
	                      0x00000000003FE000};

	uint64_t x = in;
	uint64_t t;

	x = (x | (x << 24)) & B[3];
	x = (x | (x << 12)) & B[2];
	x = (x | (x << 6)) & B[1];
	x = (x | (x << 3)) & B[0];

	for (i = 0; i < 4; i++)
	{
	#pragma HLS PIPELINE
	#pragma HLS unroll

		t = x & M[i];
		x ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);
	}

	return x & GFMASK;
}

/* input: field element in, m */
/* return: (in^2)*m */
static inline gf gf_sqmul_kernel29_2(gf in, gf m)
{
	int i;

	uint64_t x;
	uint64_t t0;
	uint64_t t1;
	uint64_t t;

	const uint64_t M[] = {0x0000001FF0000000,
	                      0x000000000FF80000,
	                      0x000000000007E000};

	t0 = in;
	t1 = m;

	x = (t1 << 6) * (t0 & (1 << 6));

	t0 ^= (t0 << 7);

	x ^= (t1 * (t0 & (0x04001)));
	x ^= (t1 * (t0 & (0x08002))) << 1;
	x ^= (t1 * (t0 & (0x10004))) << 2;
	x ^= (t1 * (t0 & (0x20008))) << 3;
	x ^= (t1 * (t0 & (0x40010))) << 4;
	x ^= (t1 * (t0 & (0x80020))) << 5;

	for (i = 0; i < 3; i++)
	{
	// #pragma HLS PIPELINE
	// #pragma HLS unroll

		t = x & M[i];
		x ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);
	}

	return x & GFMASK;
}

/* input: field element in, m */
/* return: ((in^2)^2)*m */
static inline gf gf_sq2mul_kernel29_2(gf in, gf m)
{
	int i;

	uint64_t x;
	uint64_t t0;
	uint64_t t1;
	uint64_t t;

	const uint64_t M[] = {0x1FF0000000000000,
		              0x000FF80000000000,
		              0x000007FC00000000,
	                      0x00000003FE000000,
	                      0x0000000001FE0000,
	                      0x000000000001E000};

	t0 = in;
	t1 = m;

	x = (t1 << 18) * (t0 & (1 << 6));

	t0 ^= (t0 << 21);

	x ^= (t1 * (t0 & (0x010000001)));
	x ^= (t1 * (t0 & (0x020000002))) << 3;
	x ^= (t1 * (t0 & (0x040000004))) << 6;
	x ^= (t1 * (t0 & (0x080000008))) << 9;
	x ^= (t1 * (t0 & (0x100000010))) << 12;
	x ^= (t1 * (t0 & (0x200000020))) << 15;

	for (i = 0; i < 6; i++)
	{
	// #pragma HLS PIPELINE
	// #pragma HLS unroll
		t = x & M[i];
		x ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);
	}

	return x & GFMASK;
}




gf gf_frac_kernel29_2(gf den, gf num)
{
	gf tmp_11;
	gf tmp_1111;
	gf out;

	tmp_11 = gf_sqmul_kernel29_2(den, den); // ^11
	tmp_1111 = gf_sq2mul_kernel29_2(tmp_11, tmp_11); // ^1111
	out = gf_sq2_kernel29_2(tmp_1111);
	out = gf_sq2mul_kernel29_2(out, tmp_1111); // ^11111111
	out = gf_sq2_kernel29_2(out);
	out = gf_sq2mul_kernel29_2(out, tmp_1111); // ^111111111111

	return gf_sqmul_kernel29_2(out, num); // ^1111111111110 = ^-1
}


gf gf_inv_kernel29_2(gf den)
{
	return gf_frac_kernel29_2(den, ((gf) 1));
}



gf eval_inner29_2(gf *f, gf a)
{
        int i;
        gf r;

        r = f[ SYS_T ];

        for (i = SYS_T-1; i >= 0; i--)
        {
//		#pragma HLS PIPELINE II=3
//		#pragma HLS unroll factor=2
                r = gf_mul_kernel29_2(r, a) ^ f[i];
        }

        return r;
}


void synd_kernel29_2(gf *out_out, gf *f_in, gf *L_in, unsigned char *r_in)
{

	#pragma HLS INTERFACE m_axi     port=out_out  offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=f_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=L_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=r_in     offset=slave bundle=gmem1
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
	gf local_L[SYS_N];
	gf tmp_mul_1, tmp_mul_2;
	unsigned char local_r[MAT_COLS];

	gf e_mat[SYS_N];

	// #pragma HLS ARRAY_PARTITION variable=local_out cyclic factor=4 //64
	// #pragma HLS ARRAY_PARTITION variable=local_L cyclic factor=4 //8
	// #pragma HLS ARRAY_PARTITION variable=e_mat cyclic factor=4 //8
	// #pragma HLS ARRAY_PARTITION variable=local_f cyclic factor=4 //8

	//READ into local vars

	LOOP_LOAD_FROM_BRAM_F:for (uint i=0;i<=SYS_T;i++){
	#pragma HLS PIPELINE II=1
		local_f[i] = *(f_in+i);
	}

	LOOP_LOAD_FROM_BRAM_L:for (uint i=1*SYS_N/29;i<2*SYS_N/29;i++){
	#pragma HLS PIPELINE II=1
	// #pragma HLS unroll factor=4
		local_L[i] = *(L_in+i);
	}

	LOOP_LOAD_FROM_BRAM_R:for (uint i=1*MAT_COLS/29;i<2*MAT_COLS/29;i++){
	#pragma HLS PIPELINE II=1
	// #pragma HLS unroll factor=2
		local_r[i] = *(r_in+i);
	}

	//READ into local vars END
	LOOP_EVAL:
	for(uint i=1*SYS_N/29; i <2*SYS_N/29; i++){
//	#pragma HLS PIPELINE
//	#pragma HLS unroll factor=2
		e_mat[i] = eval_inner29_2(local_f, local_L[i]);
	}



	LOOP_MAIN_OUTER:
	for (uint i = 1*SYS_N/29; i < 2*SYS_N/29; i++) //12
	{
//	#pragma HLS pipeline

		c = (local_r[i>>3] >> (i%8)) & 1;
		e_inv = gf_inv_kernel29_2(gf_mul_kernel29_2(e_mat[i],e_mat[i]));

		LOOP_MAIN_INNER:
		for (uint j = 0; j < 2*SYS_T; j++)//8
		{
//		#pragma HLS DEPENDENCE inter variable=local_out false
//		#pragma HLS PIPELINE
//		#pragma HLS unroll factor=2

			if(i==1*SYS_N/29){
				local_out[j] = gf_mul_kernel29_2(e_inv, c);
			}else{
				local_out[j] ^= gf_mul_kernel29_2(e_inv, c);

			}
			e_inv = gf_mul_kernel29_2(e_inv, local_L[i]);

		}
	}

	//WRITE to local memory

	LOOP_WRITE_TO_BRAM_OUT:for(uint i=0;i<(2*SYS_T);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		*(out_out+i) = local_out[i];
	}
}
