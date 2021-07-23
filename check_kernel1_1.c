#include "params.h"
#include "gf.h"
#include "crypto_kem.h"
#include <stdlib.h>
#include <string.h>

void check_kernel1_1(uint16_t *check_out, int *wcheck_in, gf *scheck_in, gf *s_cmpcheck_in)
{

	#pragma HLS INTERFACE m_axi     port=check_out       offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=wcheck_in       offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=scehck_in       offset=slave bundle=gmem2
	#pragma HLS INTERFACE m_axi     port=s_cmpcheck_in   offset=slave bundle=gmem3
	#pragma HLS INTERFACE s_axilite port=check_out            bundle=control
	#pragma HLS INTERFACE s_axilite port=wcheck_in            bundle=control
	#pragma HLS INTERFACE s_axilite port=scheck_in            bundle=control
	#pragma HLS INTERFACE s_axilite port=s_cmpcheck_in        bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		      bundle=control


	gf local_scheck[SYS_T*2];
	gf local_s_cmpcheck[SYS_T*2];
	uint16_t local_check;
	int local_wcheck;


	//Read the w input value from the host side
	local_wcheck = *wcheck_in;

	LOOP_LOAD_FROM_BRAM_SCHECK:for (int i=0;i<(SYS_T*2);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_scheck[i] = *(scheck_in+i);
	}

	//TODO maybe combine and read both at the same loop?
	LOOP_LOAD_FROM_BRAM_S_CMPCHECK:for (int i=0;i<(SYS_T*2);i++){
	#pragma HLS PIPELINE II=1
	#pragma HLS unroll factor=2
		local_s_cmpcheck[i] = *(s_cmpcheck_in+i);
	}


	local_check = local_wcheck;
	local_check ^= SYS_T;

	for (int i = 0; i < SYS_T*2; i++){
		#pragma HLS PIPELINE II=1
		local_check |= local_scheck[i] ^ local_s_cmpcheck[i];
	}

	local_check -= 1;
	local_check >>= 15;

	local_check ^= 1;

	//Pass he output vlue to the host side
	*check_out = local_check;




}
