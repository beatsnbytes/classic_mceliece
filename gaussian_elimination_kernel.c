
#include "pk_gen.h"
#include "params.h"
#include <stdlib.h>
#include <string.h>

//Using the same pointer for taking input and returning results

void gaussian_elimination_kernel(unsigned char *mat_in, unsigned char *mat_out, unsigned int *fail_flag)
{
	#pragma HLS INTERFACE m_axi     port=mat_in     offset=slave bundle=gmem0
	#pragma HLS INTERFACE m_axi     port=mat_out    offset=slave bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=fail_flag  offset=slave bundle=gmem0
    #pragma HLS INTERFACE s_axilite port=mat_in               bundle=control
	#pragma HLS INTERFACE s_axilite port=mat_out              bundle=control
	#pragma HLS INTERFACE s_axilite port=fail_flag            bundle=control
	#pragma HLS INTERFACE s_axilite port=return 		      bundle=control


//	unsigned int i, k, row;
//	unsigned int j, c;
//
	int i, k, row;
	int j, c;

	unsigned char mask, tmpVal, cmpVal;

	unsigned char tmpRow[MAT_COLS];
    unsigned char localMat[MAT_ROWS][MAT_COLS]; // Local memory to store input matrix

	#pragma HLS ARRAY_PARTITION variable=localMat cyclic factor=128 dim=2
	#pragma HLS ARRAY_PARTITION variable=tmpRow cyclic factor=128 //76

//	#pragma HLS RESOURCE variable=tmpRow core=RAM_2P_LUTRAM
//	#pragma HLS RESOURCE variable=localMat core=RAM_2P_LUTRAM




	LOOP_READ_FROM_DRAM_PK:
	for(i=0;i<MAT_ROWS;i++){
		for(j=0;j<MAT_COLS;j++){
			#pragma HLS PIPELINE II=1
			#pragma HLS unroll factor=4
			localMat[i][j] = *(mat_in+i*MAT_COLS+j);
		}
	}


	*fail_flag=0;
	OUTER_LOOP_1:for (i = 0; i < MAT_ROWS/8; i++)
	{
	OUTER_LOOP_2:for (j = 0; j < 8; j++)
		{
			row = (i<<3) + j;
			TMP_ROW_CONSTRUCTION_LOOP2:for(c=0;c<MAT_COLS;c++){
				#pragma HLS dependence variable=tmpRow inter false
				#pragma HLS unroll factor=128
				#pragma HLS PIPELINE II=1

				if (row>=1){
					localMat[row-1][c] = tmpRow[c]; //Store the tmpRow from previous iteration
				}
				tmpRow[c] = localMat[row][c]; //Load the tmpRow from current iteration
			}

			OUTER_LOOP_FWD_ELIM:for (k = row+1; k < MAT_ROWS; k++)
			{
				#pragma HLS dependence variable=tmpRow inter false
				#pragma HLS LOOP_TRIPCOUNT min=1 max=1643
//				#pragma HLS unroll factor=2
//				#pragma HLS PIPELINE

				mask = tmpRow[i] ^ localMat[k][i];
				mask >>= j;
				mask &= 1;
				mask = -mask;

				INNER_LOOP_FWD_ELIM:for (c = 0; c < MAT_COLS; c++)
				{
				#pragma HLS dependence variable=tmpRow inter false
				#pragma HLS PIPELINE
				#pragma HLS unroll factor=128
					tmpRow[c] ^= localMat[k][c] & mask;
				}
			}

			if ( (( tmpRow[i] >> j) & 1) == 0 ) // return if not systematic
			{
				*fail_flag = 1;
				return;
			} else {
				OUTER_LOOP_BACK_SUB:for (k = 0; k < MAT_ROWS; k++)
				{
//				#pragma HLS dependence variable=localMat inter false//not before stay as is
				#pragma HLS unroll factor=4
				#pragma HLS PIPELINE
					if (k != row)
					{
						mask = localMat[k][i] >> j;
					    (mask&1)==1?(mask=255):(mask=0);

						INNER_LOOP_BACK_SUB:for (c = 0; c < MAT_COLS; c++){
//						#pragma HLS dependence variable=tmpRow inter false
//						#pragma HLS dependence variable=localMat inter false
//						#pragma HLS PIPELINE
//						#pragma HLS unroll factor=8

							localMat[k][c] ^= tmpRow[c] & mask;
						}
					}
				}
			}
		}
	}


	LOOP_WRITE_TO_DRAM_ALT:
	for(unsigned int i=0;i<MAT_ROWS;i++){
		for(unsigned int j=0;j<MAT_COLS;j++){
		#pragma HLS PIPELINE II=1
		#pragma HLS unroll factor=2
			*(mat_out + i*MAT_COLS + j) = localMat[i][j];
		}
	}


}
