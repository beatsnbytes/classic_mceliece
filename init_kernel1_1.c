#include "params.h"
#include "gf.h"
#include "crypto_kem.h"
#include <stdlib.h>
#include <string.h>


uint16_t load_gf_init_kernel(const unsigned char *src)
{
	uint16_t a;

	a = src[1];
	a <<= 8;
	a |= src[0];

	return a & GFMASK;
}

void init_kernel1_1(gf *ginit_out, unsigned char *rinit_out, const unsigned char *cinit_in, unsigned char *skinit_in)
{


	#pragma HLS INTERFACE m_axi     port=ginit_out   	offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=rinit_out   	offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=cinit_in   	offset=slave bundle=gmem2
	#pragma HLS INTERFACE m_axi     port=skinit_in      offset=slave bundle=gmem3
	#pragma HLS INTERFACE s_axilite port=ginit_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=rinit_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=cinit_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=skinit_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		         bundle=control

	#pragma HLS inline recursive


	//TODO optimize and partition arrays
	unsigned char local_rinit[ SYS_N/8 ];
	gf local_ginit[ SYS_T+1 ];
	unsigned char local_cinit[SYND_BYTES];
	//TODO maybe needing untill 2*SYS_T....try on functional correctness
	unsigned char local_skinit[crypto_kem_SECRETKEYBYTES];
	int j=0;

	#pragma HLS ARRAY_PARTITION variable=local_skinit cyclic factor=16
	#pragma HLS ARRAY_PARTITION variable=local_ginit cyclic factor=15
//	#pragma HLS ARRAY_PARTITION variable=local_rinit cyclic factor=2


	LOOP_LOAD_FROM_BRAM_CINIT:for (int i=0;i<(SYND_BYTES);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
//		local_cinit[i] = *(cinit_in+i);
		local_rinit[i] = *(cinit_in+i);
	}

	LOOP_LOAD_FROM_BRAM_SKINIT:for (int i=0;i<crypto_kem_SECRETKEYBYTES;i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_skinit[i] = *(skinit_in+i);
	}

//	LOOP_INIT_R_1:
//	for (int i = 0; i < SYND_BYTES; i++) {
//	#pragma HLS PIPELINE II=1
//		local_rinit[i] = local_cinit[i];
//	}

	LOOP_INIT_R_2:
	for (int i = SYND_BYTES; i < SYS_N/8; i++) {
	#pragma HLS PIPELINE
		local_rinit[i] = 0;
	}

	LOOP_MAIN_COMP:
	for (int i = 0; i < SYS_T; i++) {
	#pragma HLS PIPELINE
		local_ginit[i] = load_gf_init_kernel(&local_skinit[j]);
		j += 2;
	}
	local_ginit[ SYS_T ] = 1;



	LOOP_WRITE_TO_BRAM_GOUT:for(int i=0;i<=(SYS_T);i++){
	#pragma HLS PIPELINE II=1
		*(ginit_out+i) = local_ginit[i];
	}

	LOOP_WRITE_TO_BRAM_ROUT:for(int i=0;i<(SYS_N/8);i++){
	#pragma HLS PIPELINE II=1
		*(rinit_out+i) = local_rinit[i];
	}


}
