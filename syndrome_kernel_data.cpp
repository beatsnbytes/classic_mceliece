#include "params.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ap_int.h>
#include "hls_stream.h"

void read_function_pk(unsigned char * pk_in,  hls::stream<unsigned char> &stream_pk)
{

//	LOOP_LOAD_FROM_BRAM_PK:
//	for(int i=0;i<MAT_ROWS;i++){
//		for(int j=0;j<PK_ROW_BYTES;j++){
//			#pragma HLS PIPELINE II=1
//			stream_pk.write(*(pk_in+i*PK_ROW_BYTES+j));
//		}
//	}

		LOOP_READ_FROM_DRAM_ALT:
		    for (unsigned int itr = 0, i = 0, j = 0; itr < MAT_ROWS * PK_ROW_BYTES; itr++, j++) {
		       #pragma HLS LOOP_TRIPCOUNT min=768*340 max=768*340
		       #pragma HLS PIPELINE II=1
		        if (j == PK_ROW_BYTES) {
		            j = 0;
		            i++;
		        }
		        	stream_pk.write(*(pk_in+i*PK_ROW_BYTES+j));
		    }



}


void read_function_e(unsigned char * e_in, hls::stream<unsigned char> &stream_e)
{


	LOOP_LOAD_FROM_BRAM_E:for(unsigned int i=0;i<MAT_COLS;i++){
		#pragma HLS PIPELINE II=1
		stream_e.write(*(e_in+i));
	}

}


void compute_function(hls::stream<unsigned char> &stream_pk, hls::stream<unsigned char> &stream_e, hls::stream<unsigned char> &stream_s){

	unsigned char local_s[SYND_BYTES];
	unsigned char b, row[MAT_COLS];
	unsigned char local_pk[MAT_ROWS][PK_ROW_BYTES];
	unsigned char local_e[MAT_COLS];

	#pragma HLS ARRAY_PARTITION variable=row cyclic factor=16
	#pragma HLS ARRAY_PARTITION variable=local_e cyclic factor=16
	#pragma HLS ARRAY_PARTITION variable=local_s cyclic factor=24
	#pragma HLS ARRAY_PARTITION variable=local_pk cyclic factor=16 dim=2


	LOOP_INIT_S:for (unsigned int i = 0; i < SYND_BYTES; i++){
		#pragma HLS PIPELINE
		#pragma HLS unroll factor=24
		local_s[i] = 0;
	}


	LOOP_MAIN:
	for (int i = 0; i < PK_NROWS; i++)
	{
//	#pragma HLS DEPENDENCE variable=row inter RAW true
//	#pragma HLS PIPELINE



		 LOOP_ROW_MAT:
		 for ( uint j = 0; j <(MAT_COLS); j++) {
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=16

			 if(j<(MAT_COLS - PK_ROW_BYTES)){
				 row[j] = 0;
			 }else{

//				 row[j] = local_pk[i][j-(MAT_COLS - PK_ROW_BYTES)];

				 row[j] = stream_pk.read();
			 }

		 }


		row[i>>3] |= 1 << (i%8);

		b = 0;
		LOOP_B_COMPUTE:for (uint j = 0; j < MAT_COLS; j++){
			#pragma HLS PIPELINE
			#pragma HLS unroll factor=16
			//read the local_e variables from the stream only for the first iter
			if (i==0){
				local_e[j] = stream_e.read();
			}
			b ^= row[j] & local_e[j];
		}

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;



		local_s[ i>>3 ] |= (b << (i%8));

//		Write result to output stream
		if(i%8==7){
			stream_s.write(local_s[(i>>3)]);
		}

	}

}




void write_function(unsigned char *s_out, hls::stream<unsigned char> &stream_s){


	LOOP_WRITE_TO_BRAM_R:for (unsigned int i=0;i<SYND_BYTES;i++){
		#pragma HLS PIPELINE
//		#pragma HLS unroll factor=4
		*(s_out+i) = stream_s.read();
	}
}

void syndrome_kernel_data(unsigned char *pk_in, unsigned char *e_in, unsigned char *s_out)
{
	#pragma HLS INTERFACE m_axi     port=pk_in  offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi     port=e_in   offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=s_out  offset=slave bundle=gmem2
    #pragma HLS INTERFACE s_axilite port=pk_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=e_in                bundle=control
	#pragma HLS INTERFACE s_axilite port=s_out               bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		     bundle=control

    static hls::stream<unsigned char> stream_pk;
    static hls::stream<unsigned char> stream_e;
    static hls::stream<unsigned char> stream_s;

	#pragma HLS STREAM variable = stream_pk depth = 32
	#pragma HLS STREAM variable = stream_e depth = 32
	#pragma HLS STREAM variable = stream_s depth = 32


	#pragma HLS DATAFLOW


	read_function_pk(pk_in, stream_pk);
	read_function_e(e_in, stream_e);
	compute_function(stream_pk, stream_e, stream_s);
	write_function(s_out, stream_s);


}
