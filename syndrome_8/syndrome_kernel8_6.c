#include "../params.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//#include "ap_cint.h"

void syndrome_kernel8_6(unsigned char *pk_in, unsigned char *e_in, unsigned char *s_out)
{
	#pragma HLS INTERFACE m_axi     port=pk_in    offset=slave bundle=gmem7
	#pragma HLS INTERFACE m_axi     port=e_in     offset=slave bundle=gmem8
	#pragma HLS INTERFACE m_axi     port=s_out    offset=slave bundle=gmem6
    #pragma HLS INTERFACE s_axilite port=pk_in                 bundle=control
	#pragma HLS INTERFACE s_axilite port=e_in                  bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out                 bundle=control
	#pragma HLS INTERFACE s_axilite port=return    		       bundle=control

	unsigned char b, row[MAT_COLS];


	unsigned char local_pk[MAT_ROWS/8][PK_ROW_BYTES];
	//todo canbe synd/2
	unsigned char local_s[SYND_BYTES/8];
	unsigned char local_e[MAT_COLS];

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=64
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=64
	#pragma HLS ARRAY_PARTITION variable=local_s cyclic factor=12
	#pragma HLS ARRAY_PARTITION variable=local_pk cyclic factor=64 dim=2



	LOOP_LOAD_FROM_BRAM_PK:
	for(int i=0;i<MAT_ROWS/8;i++){
		for(int j=0;j<PK_ROW_BYTES;j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=4
			local_pk[i][j] = *(pk_in+(5*PK_ROW_BYTES*MAT_ROWS/8)+i*PK_ROW_BYTES+j);
//			local_pk[i][j] = *(pk_in+(MAT_COLS*MAT_ROWS/2)+i*MAT_COLS+j + PK_NROWS/8);
		}
	}


	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS;i++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=2
		local_e[i] = *(e_in+i);
	}

		LOOP_INIT_S:for (unsigned int i = 0; i < SYND_BYTES/8; i++){
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=12
//			local_s[i+SYND_BYTES/4] = 0;
			local_s[i] = 0;
		}



	LOOP_MAIN:
	for (int i = 0; i < PK_NROWS/8; i++)
	{
	#pragma HLS PIPELINE

		 LOOP_ROW_MAT:
		 for ( uint j = 0; j <(MAT_COLS); j++) {

			#pragma HLS PIPELINE
			#pragma HLS unroll factor=64

			 if(j<(MAT_COLS - PK_ROW_BYTES)){
				 row[j] = 0;
			 }else{

				 row[j] = local_pk[i][j-(MAT_COLS - PK_ROW_BYTES)];
			 }

		 }


		row[(i>>3)+5*SYND_BYTES/8] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=64

			b ^= row[j] & local_e[j];
		}

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

//		local_s[(i>>3)+SYND_BYTES/4] |= (b << (i%8));
		local_s[(i>>3)] |= (b << (i%8));

	}


	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=5*SYND_BYTES/8;i<6*SYND_BYTES/8;i++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=4
//		*(s_out+i) = local_s[i];
		*(s_out+i) = local_s[i-5*SYND_BYTES/8];
	}

}
