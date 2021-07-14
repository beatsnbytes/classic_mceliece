#include "../params.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hls_stream.h"
#include "ap_int.h"


void read_function_pk2_2(hls::stream<data_packed_pk> &stream_pk, data_packed_pk * pk_in)
{


	LOOP_LOAD_FROM_BRAM_PK:
	for(int i=((MAT_ROWS/2)*PK_ROW_BYTES/PACK_FACTOR_PK); i<((MAT_ROWS)*PK_ROW_BYTES/PACK_FACTOR_PK); i++){
		stream_pk << *(pk_in+i);
	}

}


void read_function_e2_2(hls::stream<data_packed_e> &stream_e, data_packed_e * e_in)
{
	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS/PACK_FACTOR_E;i++){
		#pragma HLS PIPELINE II=1
		stream_e << *(e_in+i);

	}

}



void compute_function2_2(hls::stream<unsigned char> &stream_s, hls::stream<data_packed_pk> &stream_pk, hls::stream<data_packed_e> &stream_e){

	unsigned char local_s[SYND_BYTES/2][9];
	unsigned char b, row[MAT_COLS+PACK_FACTOR_PK];
	unsigned char local_e[MAT_COLS];
	data_packed_pk packed_pk_bundle;
	data_packed_e packed_e;

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=62
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=62
	#pragma HLS ARRAY_PARTITION variable=local_s cyclic factor=9 dim=2


	LOOP_INIT_S:for (uint i = 0; i <SYND_BYTES/2; i++){
		#pragma HLS PIPELINE II=1
		local_s[i][0] = 0;
	}


	LOOP_READ_FROM_STREAM_E:
	for (int i=0;i<MAT_COLS; i=i+PACK_FACTOR_E){
		#pragma HLS DEPENDENCE variable=local_e inter false
		#pragma HLS PIPELINE II=1
		stream_e >> packed_e;

		LOOP_INNER_UNPACK_E:
		for(int k=0; k<PACK_FACTOR_E; k++){
			local_e[i+k] = packed_e.packed_values[k];
		}


	}

	int done=1;
	int prev_tail=0;
	int cnt=0;

	LOOP_MAIN:
	for (int i = 0; i < PK_NROWS/2; i++)
	{
//	#pragma HLS DEPENDENCE variable=row inter RAW true
//	#pragma HLS DEPENDENCE variable=local_s inter RAW true
//	#pragma HLS PIPELINE


		LOOP_INIT_ROW:
		for(int j=0; j<(MAT_COLS - PK_ROW_BYTES); j++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=62
			row[j]=0;
		}

		done=0;
		cnt=prev_tail;
		 LOOP_ROW_MAT:
		 for ( int j = (MAT_COLS - PK_ROW_BYTES); j <MAT_COLS; j=j+PACK_FACTOR_PK) {
			#pragma HLS PIPELINE
			 if(j<MAT_COLS-prev_tail){
				 if(prev_tail!=0 && done==0){
					 LOOP_TAIL:
//					 		 for(int l=0; l<prev_tail; l++){
					 for(int l=0; l<PACK_FACTOR_PK; l++){
						 #pragma HLS PIPELINE
						 if (l<prev_tail){
							 row[MAT_COLS - PK_ROW_BYTES+l] = row[MAT_COLS+l];
						 } else {
							 break;
						 }

					 }
					 done=1;
				 }
				 //read a new bundle
				 stream_pk >> packed_pk_bundle;

				 LOOP_INNER_UNPACK_PK:
				 for(int k=0; k<PACK_FACTOR_PK; k++){
					#pragma HLS PIPELINE
					 row[j+prev_tail+k] = packed_pk_bundle.packed_values[k]; //+prev_tail
//					 printf("\n index = %d\n", curr_idx);
					 cnt++;
				 }
			 }
		 }


		 prev_tail= cnt-PK_ROW_BYTES;



		row[(i>>3)+SYND_BYTES/2] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=62
			b ^= row[j] & local_e[j];
		}

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;


//try with local variables to decouple different iterations
//		local_s[ i>>3 ] |= (b << (i%8));
		local_s[ i>>3 ][(i%8)+1] = local_s[ i>>3 ][(i%8)] | (b << (i%8)); //the end result is always at the local_s[n][7] element

	}


	LOOP_WRITE_STREAM:for (unsigned int i=0;i<SYND_BYTES/2;i++){
		#pragma HLS PIPELINE II=1
		stream_s << local_s[i][8];

	}

}



void write_function2_2(unsigned char *s_out, hls::stream<unsigned char> &stream_s){


	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=SYND_BYTES/2;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE II=1
		stream_s >> *(s_out+i);

	}
}

void syndrome_kernel_dataflow2_2(data_packed_pk *pk_in, data_packed_e *e_in, unsigned char *s_out)
{
	#pragma HLS DATAFLOW

	#pragma HLS INTERFACE m_axi bundle=gmem3 port=pk_in offset=slave//num_read_outstanding=3000
	#pragma HLS INTERFACE m_axi port=e_in   offset=slave bundle=gmem4
	#pragma HLS INTERFACE m_axi port=s_out  offset=slave bundle=gmem5
    #pragma HLS INTERFACE s_axilite bundle=control port=pk_in
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control


    static hls::stream<data_packed_pk> stream_pk;
    static hls::stream<data_packed_e> stream_e;
    static hls::stream<unsigned char> stream_s;

	#pragma HLS STREAM variable=stream_pk depth=2
	#pragma HLS STREAM variable=stream_e depth=2
	#pragma HLS STREAM variable=stream_s depth=2


	read_function_pk2_2(stream_pk, pk_in);
	read_function_e2_2(stream_e, e_in);
	compute_function2_2(stream_s, stream_pk, stream_e);
	write_function2_2(s_out, stream_s);

}
