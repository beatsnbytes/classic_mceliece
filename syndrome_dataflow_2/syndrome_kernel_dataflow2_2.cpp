#include "../params.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hls_stream.h"
#include "ap_int.h"


void read_function_pk2_2(hls::stream<ap_uint<PACK_BITWIDTH_PK>> &stream_pk, ap_uint<PACK_BITWIDTH_PK> * pk_in)
{


	LOOP_LOAD_FROM_BRAM_PK:
	for(uint i=0;i<MAT_ROWS;i++){
		for(uint j=0;j<PK_ROW_BYTES/PACK_FACTOR_PK;j++){
			#pragma HLS PIPELINE II=1
			stream_pk << *(pk_in+i*PK_ROW_BYTES/PACK_FACTOR_PK+j);
		}
	}

}


void read_function_e2_2(hls::stream<ap_uint<PACK_BITWIDTH_E>> &stream_e, ap_uint<PACK_BITWIDTH_E> * e_in)
{
	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS/PACK_FACTOR_E;i++){
		#pragma HLS PIPELINE II=1
		stream_e << *(e_in+i);

	}

}



void compute_function2_2(hls::stream<unsigned char> &stream_s, hls::stream<ap_uint<PACK_BITWIDTH_PK>> &stream_pk, hls::stream<ap_uint<PACK_BITWIDTH_E>> &stream_e){

	unsigned char local_s[SYND_BYTES];
	unsigned char b, row[MAT_COLS];
	unsigned char local_e[MAT_COLS];

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=64
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=64


	LOOP_INIT_S:for (uint i = 0; i <SYND_BYTES; i++){
		#pragma HLS PIPELINE II=1
		local_s[i] = 0;
	}


	ap_uint<PACK_BITWIDTH_E> packed_e;
	LOOP_READ_FROM_STREAM_E:
	for (int i=0;i<MAT_COLS; i=i+PACK_FACTOR_E){
		#pragma HLS DEPENDENCE variable=local_e inter false
		#pragma HLS PIPELINE II=1
		stream_e >> packed_e;

		LOOP_INNER_UNPACK_E:
		for(int k=0; k<PACK_FACTOR_E; k++){
			local_e[i+k] = (packed_e >> (k<<3) ) & 0xff;
		}


	}


	LOOP_MAIN:
	for (int i = 0; i < PK_NROWS; i++)
	{
//	#pragma HLS DEPENDENCE variable=row inter RAW true
	#pragma HLS PIPELINE


		LOOP_INIT_ROW:
		for(int j=0; j<(MAT_COLS - PK_ROW_BYTES); j++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=64
			row[j]=0;
		}


		ap_uint<PACK_BITWIDTH_PK> ext32b_row;

		 LOOP_ROW_MAT:
		 for ( uint j = (MAT_COLS - PK_ROW_BYTES); j <MAT_COLS; j=j+PACK_FACTOR_PK) {
			#pragma HLS DEPENDENCE variable=row inter false
			#pragma HLS DEPENDENCE variable=stream_pk inter false
			#pragma HLS PIPELINE II=1
				 //The stream provides 32b values
				 stream_pk >> ext32b_row;

				 LOOP_INNER_UNPACK_PK:
				 for(int k=0; k<PACK_FACTOR_PK; k++){
					 row[j+k] = (ext32b_row >> (k<<3) ) & 0xff;
				 }

			 }



		row[i>>3] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=64
			b ^= row[j] & local_e[j];
		}

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;



		local_s[ i>>3 ] |= (b << (i%8));

	}


	LOOP_WRITE_STREAM:for (unsigned int i=0;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE II=1
		stream_s << local_s[i];

	}

}




void write_function2_2(unsigned char *s_out, hls::stream<unsigned char> &stream_s){


	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=0;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE
		stream_s >> *(s_out+i);

	}
}

void syndrome_kernel_dataflow2_2(ap_uint<PACK_BITWIDTH_PK> *pk_in, ap_uint<PACK_BITWIDTH_E> *e_in, unsigned char *s_out)
{
	#pragma HLS DATAFLOW

	#pragma HLS INTERFACE m_axi depth=10000 bundle=gmem port=pk_in offset=slave//depth=4080 num_read_outstanding=3000
	#pragma HLS INTERFACE m_axi port=e_in   offset=slave bundle=gmem1 depth=436//depth=109
	#pragma HLS INTERFACE m_axi port=s_out  offset=slave bundle=gmem2 depth=96//depth=3
    #pragma HLS INTERFACE s_axilite bundle=control port=pk_in
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control


    static hls::stream<ap_uint<PACK_BITWIDTH_PK>> stream_pk;
    static hls::stream<ap_uint<PACK_BITWIDTH_E>> stream_e;
    static hls::stream<unsigned char> stream_s;

	#pragma HLS STREAM variable=stream_pk depth=2
	#pragma HLS STREAM variable=stream_e depth=2
	#pragma HLS STREAM variable=stream_s depth=2


	read_function_pk2_2(stream_pk, pk_in);
	read_function_e2_2(stream_e, e_in);
	compute_function2_2(stream_s, stream_pk, stream_e);
	write_function2_2(s_out, stream_s);

}
