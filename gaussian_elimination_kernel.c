#include "pk_gen.h"
#include "params.h"
#include <stdlib.h>
#include <string.h>
//#include "ap_cint.h"

//Using the same pointer for taking input and returning results

void gaussian_elimination_kernel(unsigned char *mat_in, unsigned char *mat_out)
{
	#pragma HLS INTERFACE m_axi     port=mat_in   offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi     port=mat_out  offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=mat_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=mat_out              bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		      bundle=control


	unsigned int i, k, row;
	unsigned int j, c;

	unsigned char mask, tmpVal, cmpVal;

	unsigned char tmpRow[MAT_COLS];
    unsigned char localMat[MAT_ROWS][MAT_COLS]; // Local memory to store input matrix

	#pragma HLS ARRAY_PARTITION variable=localMat cyclic factor=218 dim=2
	#pragma HLS ARRAY_PARTITION variable=tmpRow cyclic factor=218


//	LOOP_READ_FROM_DRAM_PK:
//	for(i=0;i<MAT_ROWS;i++){
//		for(j=0;j<MAT_COLS;j++){
////		#pragma HLS dependence variable=localMat inter false
//		#pragma HLS PIPELINE II=1
//		#pragma HLS unroll factor=2
//			localMat[i][j] = *(mat_in+i*MAT_COLS+j);
//		}
//	}

	LOOP_READ_FROM_DRAM_PK:
	for(i=0;i<MAT_ROWS;i++){
		for(j=0;j<MAT_COLS;j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=4
			localMat[i][j] = *(mat_in+i*MAT_COLS+j);
		}
	}



	OUTER_LOOP_1:for (i = 0; i < MAT_ROWS/8; i++)
	{
	OUTER_LOOP_2:for (j = 0; j < 8; j++)
		{
			row = (i<<3) + j;
			TMP_ROW_CONSTRUCTION_LOOP2:for(c=0;c<MAT_COLS;c++){
				#pragma HLS dependence variable=tmpRow inter false
//				#pragma HLS dependence variable=localMat inter false
				#pragma HLS unroll factor=218
				#pragma HLS PIPELINE II=1

				if (row>=1){
					localMat[row-1][c] = tmpRow[c]; //Store the tmpRow from previous iteration
				}
				tmpRow[c] = localMat[row][c]; //Load the tmpRow from current iteration
			}

			OUTER_LOOP_FWD_ELIM:for (k = row+1; k < MAT_ROWS; k++)
			{
//				#pragma HLS dependence variable=tmpRow inter RAW
				#pragma HLS dependence variable=tmpRow inter false
				#pragma HLS LOOP_TRIPCOUNT min=1 max=767
//				#pragma HLS unroll factor=2
//				#pragma HLS PIPELINE

				mask = tmpRow[i] ^ localMat[k][i];
				mask >>= j;
				mask &= 1;
				mask = -mask;

				INNER_LOOP_FWD_ELIM:for (c = 0; c < MAT_COLS; c++)
				{
				#pragma HLS dependence variable=tmpRow inter false
				#pragma HLS PIPELINE II
				#pragma HLS unroll factor=218
					tmpRow[c] ^= localMat[k][c] & mask;
				}
			}

			if ( (( tmpRow[i] >> j) & 1) == 0 ) // return if not systematic
			{
				*mat_out = 255;
				return;
			} else {
				OUTER_LOOP_BACK_SUB:for (k = 0; k < MAT_ROWS; k++)
				{
//				#pragma HLS dependence variable=localMat inter false
				#pragma HLS unroll factor=2
				#pragma HLS PIPELINE
					if (k != row)
					{
						mask = localMat[k][i] >> j;
					    (mask&1)==1?(mask=255):(mask=0);

						INNER_LOOP_BACK_SUB:for (c = 0; c < MAT_COLS; c++){
//						#pragma HLS dependence variable=tmpRow inter false
//						#pragma HLS dependence variable=localMat inter false
//						#pragma HLS PIPELINE II
//						#pragma HLS unroll factor=64

							localMat[k][c] ^= tmpRow[c] & mask;
						}
					}
				}
			}
		}
	}


	LOOP_WRITE_TO_DRAM_ALT:
	for(i=0;i<MAT_ROWS;i++){
		for(j=0;j<MAT_COLS;j++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=2
			*(mat_out+i*MAT_COLS+j) = localMat[i][j];
		}
	}

}
