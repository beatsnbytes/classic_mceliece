#include "params.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


void syndrome_kernel(unsigned char *pk_in, unsigned char *e_in, unsigned char *s_out)
{
	#pragma HLS INTERFACE m_axi     port=pk_in  offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi     port=e_in   offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=s_out  offset=slave bundle=gmem2
    #pragma HLS INTERFACE s_axilite port=pk_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

	unsigned char b, row[MAT_COLS]; //cols=870

	unsigned char local_pk[MAT_ROWS][PK_ROW_BYTES];//677 is prime
	unsigned char local_s[SYND_BYTES]; //synd_bytes=194
	unsigned char local_e[MAT_COLS];

	int tail = PK_NROWS % 8;

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=87
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=87
	#pragma HLS ARRAY_PARTITION variable=local_s cyclic factor=4 //194
	#pragma HLS ARRAY_PARTITION variable=local_pk cyclic factor=87 dim=2


	LOOP_LOAD_FROM_BRAM_PK:
	for(int i=0;i<MAT_ROWS;i++){
		for(int j=0;j<PK_ROW_BYTES;j++){
			#pragma HLS PIPELINE ΙΙ=1
			#pragma HLS unroll factor=4
			local_pk[i][j] = *(pk_in+i*PK_ROW_BYTES+j);
		}
	}


	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS;i++){
		#pragma HLS PIPELINE ΙΙ=1
		#pragma HLS unroll factor=2
		local_e[i] = *(e_in+i);
	}


	LOOP_INIT_S:for (unsigned int i = 0; i < SYND_BYTES; i++){
		#pragma HLS PIPELINE ΙΙ=1
		#pragma HLS unroll factor=4
		local_s[i] = 0;
	}





	LOOP_MAIN:
	for (int i = 0; i < PK_NROWS; i++)
	{
//	#pragma HLS DEPENDENCE variable=row inter RAW true
//	#pragma HLS PIPELINE

		LOOP_INIT_ROW:for (unsigned int i = 0; i < MAT_COLS- PK_ROW_BYTES; i++){
			#pragma HLS PIPELINE ΙΙ=1
			#pragma HLS unroll factor=87
			row[i] = 0;
		}


		 LOOP_ROW_MAT:
		 for (int j = (MAT_COLS - PK_ROW_BYTES); j <(MAT_COLS); j++) { //mat_cols
			#pragma HLS DEPENDENCE variable=local_pk inter false
			#pragma HLS DEPENDENCE variable=row inter false
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=87

				 row[j] = local_pk[i][j-(MAT_COLS - PK_ROW_BYTES)];
		 }


		for (int j = SYS_N/8-1; j >= SYS_N/8 - PK_ROW_BYTES; j--){
		#pragma HLS DEPENDENCE variable=row inter false
		#pragma HLS PIPELINE
		#pragma HLS unroll factor=87
			row[ j ] = (row[ j ] << tail) | (row[j-1] >> (8-tail));
		}


		row[i>>3] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=87

			b ^= row[j] & local_e[j];
		}

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

		local_s[ i>>3 ] |= (b << (i%8));

	}

	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=0;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE ΙΙ=1
		#pragma HLS unroll factor=4
		*(s_out+i) = local_s[i];
	}

}
