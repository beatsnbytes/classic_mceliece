#include "params.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>


#define min(a, b) ((a < b) ? a : b)

gf gf_mul_kernel_bm(gf in0, gf in1)
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

/* input: field element in */
/* return: in^2 */
static inline gf gf_sq_kernel_bm(gf in)
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

gf gf_inv_kernel_bm(gf in)
{
	gf tmp_11;
	gf tmp_1111;

	gf out = in;

	out = gf_sq_kernel_bm(out);
	tmp_11 = gf_mul_kernel_bm(out, in); // 11

	out = gf_sq_kernel_bm(tmp_11);
	out = gf_sq_kernel_bm(out);
	tmp_1111 = gf_mul_kernel_bm(out, tmp_11); // 1111

	out = gf_sq_kernel_bm(tmp_1111);
	out = gf_sq_kernel_bm(out);
	out = gf_sq_kernel_bm(out);
	out = gf_sq_kernel_bm(out);
	out = gf_mul_kernel_bm(out, tmp_1111); // 11111111

	out = gf_sq_kernel_bm(out);
	out = gf_sq_kernel_bm(out);
	out = gf_mul_kernel_bm(out, tmp_11); // 1111111111

	out = gf_sq_kernel_bm(out);
	out = gf_mul_kernel_bm(out, in); // 11111111111

	return gf_sq_kernel_bm(out); // 111111111110
}

/* input: field element den, num */
/* return: (num/den) */
gf gf_frac_kernel_bm(gf den, gf num)
{
	return gf_mul_kernel_bm(gf_inv_kernel_bm(den), num);
}



/* the Berlekamp-Massey algorithm */
/* input: s, sequence of field elements */
/* output: out, minimal polynomial of s */
void bm_kernel1_1(gf *outbm_out, gf *sbm_in)
{

	#pragma HLS INTERFACE m_axi     port=outbm_out   offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=sbm_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE s_axilite port=outbm_out             bundle=control
	#pragma HLS INTERFACE s_axilite port=sbm_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

	#pragma HLS inline recursive


	int i;
	uint16_t N = 0;
	uint16_t L = 0;
	uint16_t mle;
	uint16_t mne;

	gf T[ SYS_T+1 ];
	gf C[ SYS_T+1 ];
	gf B[ SYS_T+1 ];

	gf b = 1, d, f;


	gf local_sbm[SYS_T*2];
	gf local_outbm[SYS_T+1];


	//

	LOOP_LOAD_FROM_BRAM_SBM:for (int i=0;i<(SYS_T*2);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_sbm[i] = *(sbm_in+i);
	}





	for (i = 0; i < SYS_T+1; i++)
		C[i] = B[i] = 0;

	B[1] = C[0] = 1;

	//TODO optimize the loops below

	LOOP_MAIN_OUTTER:
	for (N = 0; N < 2 * SYS_T; N++)
	{
	#pragma HLS PIPELINE
		d = 0;

		LOOP_INNER_1:
		//TODO problem because not known beforehand. what is the dynamic range of these values?
		for (i = 0; i <= min(N, SYS_T); i++)
			d ^= gf_mul_kernel_bm(C[i], local_sbm[ N-i]);

		mne = d; mne -= 1;   mne >>= 15; mne -= 1;
		mle = N; mle -= 2*L; mle >>= 15; mle -= 1;
		mle &= mne;

		LOOP_INNER_2:
		for (i = 0; i <= SYS_T; i++)
		#pragma HLS PIPELINE II=1
			T[i] = C[i];

		f = gf_frac_kernel_bm(b, d);

		LOOP_INNER_3:
		for (i = 0; i <= SYS_T; i++)
		#pragma HLS PIPELINE II=1
			C[i] ^= gf_mul_kernel_bm(f, B[i]) & mne;

		L = (L & ~mle) | ((N+1-L) & mle);

		LOOP_INNER_4:
		for (i = 0; i <= SYS_T; i++)
		#pragma HLS PIPELINE II=1
			B[i] = (B[i] & ~mle) | (T[i] & mle);

		b = (b & ~mle) | (d & mle);

		LOOP_INNER_5:
		for (i = SYS_T; i >= 1; i--){
			#pragma HLS PIPELINE II=1
			B[i] = B[i-1];
		}

		B[0] = 0;
	}

	for (i = 0; i <= SYS_T; i++)
		local_outbm[i] = C[ SYS_T-i ];

	//TODO maybe keep only the above loop
	LOOP_WRITE_TO_BRAM_OUTBM:for(int i=0;i<=(SYS_T);i++){
	#pragma HLS PIPELINE II=1
		*(outbm_out+i) = local_outbm[i];
	}


}

