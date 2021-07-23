#include "params.h"
#include "gf.h"
#include "crypto_kem.h"
#include <stdlib.h>
#include <string.h>


gf gf_iszero_mid_kernel(gf a)
{
	uint32_t t = a;

	t -= 1;
	t >>= 19;

	return (gf) t;
}


void mid_kernel1_1(int *wmid_out, unsigned char *emid_out, gf *imagesmid_in)
{
	#pragma HLS INTERFACE m_axi     port=wmid_out   	    offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=emid_out   	    offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=imagesmid_in   	offset=slave bundle=gmem2
	#pragma HLS INTERFACE s_axilite port=wmid_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=emid_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=imagesmid_in           bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		        bundle=control


	#pragma HLS inline recursive

	int local_wmid;
	unsigned char local_emid[SYS_N/8][9];
	gf local_imagesmid[SYS_N];
	gf t;

	#pragma HLS ARRAY_PARTITION variable=local_emid cyclic factor=16
	#pragma HLS ARRAY_PARTITION variable=local_imagesmid cyclic factor=16


	LOOP_LOAD_FROM_BRAM_IMAGESMID_IN:for (int i=0;i<(SYS_N);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_imagesmid[i] = *(imagesmid_in+i);
	}


	LOOP_INIT_EMID:
	for (int i = 0; i < (SYS_N/8); i++){
		#pragma HLS PIPELINE
//		#pragma HLS unroll factor=2
		local_emid[i][0] = 0;
	}

	//TODO add unroll and array partition
	LOOP_MAIN_COMP:
	for (int i = 0; i < SYS_N; i++)
	{
	#pragma HLS PIPELINE
	#pragma HLS unroll factor=32

		t = gf_iszero_mid_kernel(local_imagesmid[i]) & 1;

//		local_emid[ i/8 ] |= t << (i%8);
		local_emid[ i>>3 ][(i%8)+1] = local_emid[ i>>3 ][(i%8)] | (t << (i%8));

		local_wmid += t;
	}

	*wmid_out = local_wmid;

	LOOP_WRITE_TO_BRAM_EMID_OUT:for(int i=0;i<(SYS_N/8);i++){
	#pragma HLS PIPELINE II=1
		*(emid_out+i) = local_emid[i][8];
	}

}
