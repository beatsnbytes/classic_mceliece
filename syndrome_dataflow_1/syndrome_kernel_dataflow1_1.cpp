#include "/tools/Xilinx/Vivado/2020.1/include/gmp.h"
#include "../params.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hls_stream.h"
#include "ap_int.h"

#define PACK_FACTOR_PK 10
#define PACK_BITWIDTH_PK PACK_FACTOR_PK*8
#define PACK_FACTOR_E 4
#define PACK_BITWIDTH_E PACK_FACTOR_E*8

void read_function_pk1_1(hls::stream<ap_uint<PACK_BITWIDTH_PK>> &stream_pk, ap_uint<PACK_BITWIDTH_PK> * pk_in)
{


	LOOP_LOAD_FROM_BRAM_PK:
	for(uint i=0;i<MAT_ROWS;i++){
		for(uint j=0;j<PK_ROW_BYTES/PACK_FACTOR_PK;j++){
			#pragma HLS PIPELINE II=1
			stream_pk << *(pk_in+i*PK_ROW_BYTES/PACK_FACTOR_PK+j);
		}
	}

}


void read_function_e1_1(hls::stream<ap_uint<PACK_BITWIDTH_E>> &stream_e, ap_uint<PACK_BITWIDTH_E> * e_in)
{
	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS/PACK_FACTOR_E;i++){
		#pragma HLS PIPELINE II=1
		stream_e << *(e_in+i);

	}

}



void compute_function1_1(hls::stream<unsigned char> &stream_s, hls::stream<ap_uint<PACK_BITWIDTH_PK>> &stream_pk, hls::stream<ap_uint<PACK_BITWIDTH_E>> &stream_e, hls::stream<unsigned char> &stream_pk_out){

	unsigned char local_s[SYND_BYTES];
	unsigned char b, row[MAT_COLS];
	unsigned char local_e[MAT_COLS];
	ap_uint<PACK_BITWIDTH_PK> ext32b_row; //previously was inside the loop

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=32
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=32


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
		#pragma HLS unroll factor=32
			row[j]=0;
		}



		 LOOP_ROW_MAT:
		 for ( int j = (MAT_COLS - PK_ROW_BYTES); j <MAT_COLS; j=j+PACK_FACTOR_PK) {
//			#pragma HLS DEPENDENCE variable=row inter false
//			#pragma HLS DEPENDENCE variable=stream_pk inter false
			#pragma HLS PIPELINE II=1
				 //The stream provides 32b values
				 stream_pk >> ext32b_row;

				 LOOP_INNER_UNPACK_PK:
				 for(int k=0; k<PACK_FACTOR_PK; k++){
					 row[j+k] = (ext32b_row >> (k<<3) ) & 0xff;
					 stream_pk_out << row[j+k];
				 }

			 }



		row[i>>3] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=32
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




void write_function1_1(unsigned char *s_out, hls::stream<unsigned char> &stream_s){


	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=0;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE
		stream_s >> *(s_out+i);

	}

}


void write_function_pk(unsigned char *pk_out, hls::stream<unsigned char> &stream_pk_out){

	LOOP_WRITE_PK_OUT:
	for(int i=0;i<MAT_ROWS;i++){
		for(int j=0;j<PK_ROW_BYTES;j++){
			#pragma HLS PIPELINE
			stream_pk_out >> *(pk_out+i*PK_ROW_BYTES+j);
		}
	}
}


void syndrome_kernel_dataflow1_1(ap_uint<PACK_BITWIDTH_PK> *pk_in, ap_uint<PACK_BITWIDTH_E> *e_in, unsigned char *s_out, unsigned char *pk_out)
{
	#pragma HLS DATAFLOW

	#pragma HLS INTERFACE m_axi port=pk_in  offset=slave bundle=gmem  depth=10000
	#pragma HLS INTERFACE m_axi port=e_in   offset=slave bundle=gmem1 depth=10000
	#pragma HLS INTERFACE m_axi port=s_out  offset=slave bundle=gmem2 depth=90000
	#pragma HLS INTERFACE m_axi port=pk_out  offset=slave bundle=gmem3 depth=190000



    #pragma HLS INTERFACE s_axilite port=pk_in				 bundle=control
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control

	#pragma HLS INTERFACE s_axilite port=pk_out               bundle=control

	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control


    static hls::stream<ap_uint<PACK_BITWIDTH_PK>> stream_pk;
    static hls::stream<ap_uint<PACK_BITWIDTH_E>> stream_e;
    static hls::stream<unsigned char> stream_s;

    static hls::stream<unsigned char> stream_pk_out;

	#pragma HLS STREAM variable=stream_pk depth=2
	#pragma HLS STREAM variable=stream_e depth=2
	#pragma HLS STREAM variable=stream_s depth=2

	#pragma HLS STREAM variable=stream_pk_out depth=2


	read_function_pk1_1(stream_pk, pk_in);
	read_function_e1_1(stream_e, e_in);
	compute_function1_1(stream_s, stream_pk, stream_e, stream_pk_out);
	write_function1_1(s_out, stream_s);;
	write_function_pk(pk_out, stream_pk_out);

}
