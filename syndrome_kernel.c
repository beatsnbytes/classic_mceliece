#include "params.h"
#include <stdlib.h>
#include <string.h>

void syndrome_kernel(unsigned char *pk_in, unsigned char *e_in, unsigned char *s_out)
{
	#pragma HLS INTERFACE m_axi     port=pk_in  offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=e_in   offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi     port=s_out  offset=slave bundle=gmem1
    #pragma HLS INTERFACE s_axilite port=pk_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

	unsigned char b, row[MAT_COLS];
	int i, j;
	unsigned char local_pk[MAT_ROWS][PK_ROW_BYTES];
	unsigned char local_s[SYND_BYTES];
	unsigned char local_e[MAT_COLS];

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=64
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=64
	#pragma HLS ARRAY_PARTITION variable=local_s cyclic factor=64 //prev 96
	#pragma HLS ARRAY_PARTITION variable=local_pk cyclic factor=64 dim=2
//	#pragma HLS ARRAY_PARTITION variable=local_pk cyclic factor=2 dim=1

	LOOP_LOAD_FROM_BRAM_PK:
	for(i=0;i<MAT_ROWS;i++){
		for(j=0;j<PK_ROW_BYTES;j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=4
			local_pk[i][j] = *(pk_in+i*PK_ROW_BYTES+j);
		}
	}

	LOOP_LOAD_FROM_BRAM_E:for(i=0;i<MAT_COLS;i++){
		#pragma HLS unroll factor=4
		local_e[i] = *(e_in+i);
	}


	LOOP_INIT_S:for (i = 0; i < SYND_BYTES; i++){
		#pragma HLS unroll factor=64
		local_s[i] = 0;
	}

	LOOP_MAIN:for (i = 0; i < PK_NROWS; i++)
	{
		#pragma HLS PIPELINE
		#pragma HLS unroll factor=2

		LOOP_ROW_INIT:for (j = 0; j < MAT_COLS; j++)
			#pragma HLS unroll factor=64
			row[j] = 0;

		LOOP_ROW_COMPUTE:for (j = 0; j < PK_ROW_BYTES; j++){
			#pragma HLS unroll factor=64
			row[ MAT_COLS - PK_ROW_BYTES + j ] = local_pk[i][j];
		}

		row[i/8] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (j = 0; j < MAT_COLS; j++)
			#pragma HLS unroll factor=64
			b ^= row[j] & local_e[j];

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

		local_s[ i/8 ] |= (b << (i%8));

	}

	LOOP_WRITE_TO_BRAM_R:for (i=0;i<SYND_BYTES;i++){
		#pragma HLS unroll factor=16
		*(s_out+i) = local_s[i];
	}

}
