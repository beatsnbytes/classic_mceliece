#include "params.h"
#include "gf.h"
#include "crypto_kem.h"
#include <stdlib.h>
#include <string.h>



uint32_t load4_kernel_supp(const unsigned char * in)
{
	int i;
	uint32_t ret = in[3];

	for (i = 2; i >= 0; i--)
	{
		ret <<= 8;
		ret |= in[i];
	}

	return ret;
}


void store8_kernel_supp(unsigned char *out, uint64_t in)
{
	out[0] = (in >> 0x00) & 0xFF;
	out[1] = (in >> 0x08) & 0xFF;
	out[2] = (in >> 0x10) & 0xFF;
	out[3] = (in >> 0x18) & 0xFF;
	out[4] = (in >> 0x20) & 0xFF;
	out[5] = (in >> 0x28) & 0xFF;
	out[6] = (in >> 0x30) & 0xFF;
	out[7] = (in >> 0x38) & 0xFF;
}



uint64_t load8_kernel_supp(const unsigned char * in)
{
	int i;
	uint64_t ret = in[7];

	for (i = 6; i >= 0; i--)
	{
		ret <<= 8;
		ret |= in[i];
	}

	return ret;
}


void transpose_64x64_kernel_supp(uint64_t * out, uint64_t * in)
{
	int i, j, s, d;

	uint64_t x, y;
	uint64_t masks[6][2] = {
	                        {0x5555555555555555, 0xAAAAAAAAAAAAAAAA},
	                        {0x3333333333333333, 0xCCCCCCCCCCCCCCCC},
	                        {0x0F0F0F0F0F0F0F0F, 0xF0F0F0F0F0F0F0F0},
	                        {0x00FF00FF00FF00FF, 0xFF00FF00FF00FF00},
	                        {0x0000FFFF0000FFFF, 0xFFFF0000FFFF0000},
	                        {0x00000000FFFFFFFF, 0xFFFFFFFF00000000}
	                       };

	for (i = 0; i < 64; i++)
		out[i] = in[i];

	for (d = 5; d >= 0; d--)
	{
		s = 1 << d;

		for (i = 0; i < 64; i += s*2)
		for (j = i; j < i+s; j++)
		{
			x = (out[j] & masks[d][0]) | ((out[j+s] & masks[d][0]) << s);
			y = ((out[j] & masks[d][1]) >> s) | (out[j+s] & masks[d][1]);

			out[j+0] = x;
			out[j+s] = y;
		}
	}
}

/* one layer of the benes network */
static void layer_kernel_supp(uint64_t * data, uint64_t * bits, int lgs)
{
	int i, j, s;

	uint64_t d;

	s = 1 << lgs;

	for (i = 0; i < 64; i += s*2)
	for (j = i; j < i+s; j++)
	{

		d = (data[j+0] ^ data[j+s]);
		d &= (*bits++);
		data[j+0] ^= d;
		data[j+s] ^= d;
	}
}

/* input: r, sequence of bits to be permuted */
/*        bits, condition bits of the Benes network */
/*        rev, 0 for normal application; !0 for inverse */
/* output: r, permuted bits */
void apply_benes_kernel_supp(unsigned char * r, const unsigned char * bits, int rev)
{
	int i;

	const unsigned char *cond_ptr;
	int inc, low;

	uint64_t bs[64];
	uint64_t cond[64];

	//

	for (i = 0; i < 64; i++)
	{
		bs[i] = load8_kernel_supp(r + i*8);
	}

	if (rev == 0)
	{
		inc = 256;
		cond_ptr = bits;
	}
	else
	{
		inc = -256;
		cond_ptr = bits + (2*GFBITS-2)*256;
	}

	//

	transpose_64x64_kernel_supp(bs, bs);

	for (low = 0; low <= 5; low++)
	{
		for (i = 0; i < 64; i++) cond[i] = load4_kernel_supp(cond_ptr + i*4);
		transpose_64x64_kernel_supp(cond, cond);
		layer_kernel_supp(bs, cond, low);
		cond_ptr += inc;
	}

	transpose_64x64_kernel_supp(bs, bs);

	for (low = 0; low <= 5; low++)
	{
		for (i = 0; i < 32; i++) cond[i] = load8_kernel_supp(cond_ptr + i*8);
		layer_kernel_supp(bs, cond, low);
		cond_ptr += inc;
	}
	for (low = 4; low >= 0; low--)
	{
		for (i = 0; i < 32; i++) cond[i] = load8_kernel_supp(cond_ptr + i*8);
		layer_kernel_supp(bs, cond, low);
		cond_ptr += inc;
	}

	transpose_64x64_kernel_supp(bs, bs);

	for (low = 5; low >= 0; low--)
	{
		for (i = 0; i < 64; i++) cond[i] = load4_kernel_supp(cond_ptr + i*4);
		transpose_64x64_kernel_supp(cond, cond);
		layer_kernel_supp(bs, cond, low);
		cond_ptr += inc;
	}

	transpose_64x64_kernel_supp(bs, bs);

	//

	for (i = 0; i < 64; i++)
	{
		store8_kernel_supp(r + i*8, bs[i]);
	}
}



gf bitrev_kernel_supp(gf a)
{
	a = ((a & 0x00FF) << 8) | ((a & 0xFF00) >> 8);
	a = ((a & 0x0F0F) << 4) | ((a & 0xF0F0) >> 4);
	a = ((a & 0x3333) << 2) | ((a & 0xCCCC) >> 2);
	a = ((a & 0x5555) << 1) | ((a & 0xAAAA) >> 1);

	return a >> 4;
}


void support_gen_kernel1_1(gf * ssupp_out, const unsigned char *csupp_in)
{
	#pragma HLS INTERFACE m_axi     port=ssupp_out   offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=csupp_in     offset=slave bundle=gmem1
	#pragma HLS INTERFACE s_axilite port=ssupp_out             bundle=control
	#pragma HLS INTERFACE s_axilite port=csupp_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

	#pragma HLS inline recursive


	gf a;
	int i, j;
	unsigned char L[ GFBITS ][ (1 << GFBITS)/8 ];

	unsigned char local_csupp[crypto_kem_SECRETKEYBYTES];
	gf local_ssupp[SYS_N];


	LOOP_LOAD_FROM_BRAM_CSUPP:for (int i=0;i<crypto_kem_SECRETKEYBYTES;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_csupp[i] = *(csupp_in+i);
	}



	for (i = 0; i < GFBITS; i++)
		for (j = 0; j < (1 << GFBITS)/8; j++)
			L[i][j] = 0;

	for (i = 0; i < (1 << GFBITS); i++)
	{
	#pragma HLS PIPELINE
		a = bitrev_kernel_supp((gf) i);

		for (j = 0; j < GFBITS; j++)
			L[j][ i/8 ] |= ((a >> j) & 1) << (i%8);
	}

	//TODO is the local_csupp right here?
	for (j = 0; j < GFBITS; j++)
//	#pragma HLS PIPELINE
		apply_benes_kernel_supp(L[j], local_csupp, 0);

	for (i = 0; i < SYS_N; i++)
	{
//	#pragma HLS PIPELINE
		local_ssupp[i] = 0;
		for (j = GFBITS-1; j >= 0; j--)
		{
			local_ssupp[i] <<= 1;
			local_ssupp[i] |= (L[j][i/8] >> (i%8)) & 1;
		}
	}


	LOOP_WRITE_TO_BRAM_OUTBM:for(int i=0;i<SYS_N;i++){
	#pragma HLS PIPELINE II=1
		*(ssupp_out+i) = local_ssupp[i];
	}


}
