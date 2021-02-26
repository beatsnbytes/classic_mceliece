/*
   PQCgenKAT_kem.c
   Created by Bassham, Lawrence E (Fed) on 8/29/17.
   Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
   + mods from djb: see KATNOTES
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"
#include "crypto_kem.h"
#include "pk_gen.h"
#include "root.h"
#include "encrypt.h"
#include "decrypt.h"
#include "operations.h"

#include <sys/time.h>
#include <valgrind/callgrind.h>

#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_CRYPTO_FAILURE  -4

void	fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);

unsigned char entropy_input[48];
unsigned char seed[KATNUM][48];

double sum_keygen, sum_enc, sum_dec;
int times_keygen, times_enc, times_dec;

int
main()
{
    FILE                *fp_req, *fp_rsp;
    int                 ret_val;
    int i;
    unsigned char *ct = 0;
    unsigned char *ss = 0;
    unsigned char *ss1 = 0;
    unsigned char *pk = 0;
    unsigned char *sk = 0;

    struct timeval start_keygen, end_keygen, start_enc, end_enc, start_dec, end_dec;

    for (i=0; i<48; i++)
        entropy_input[i] = i;
    randombytes_init(entropy_input, NULL, 256);

    for (i=0; i<KATNUM; i++)
        randombytes(seed[i], 48);

    fp_req = fopen("kat_kem.req", "w");
    if (!fp_req)
        return KAT_FILE_OPEN_ERROR;

    for (i=0; i<KATNUM; i++) {
        fprintf(fp_req, "count = %d\n", i);
        fprintBstr(fp_req, "seed = ", seed[i], 48);
        fprintf(fp_req, "pk =\n");
        fprintf(fp_req, "sk =\n");
        fprintf(fp_req, "ct =\n");
        fprintf(fp_req, "ss =\n\n");
    }

    fp_rsp = fopen("kat_kem.rsp", "w");
    if (!fp_rsp)
        return KAT_FILE_OPEN_ERROR;

    fprintf(fp_rsp, "# kem/%s\n\n", crypto_kem_PRIMITIVE);

    for (i=0; i<KATNUM; i++) {
        if (!ct) ct = malloc(crypto_kem_CIPHERTEXTBYTES);
        if (!ct) abort();
        if (!ss) ss = malloc(crypto_kem_BYTES);
        if (!ss) abort();
        if (!ss1) ss1 = malloc(crypto_kem_BYTES);
        if (!ss1) abort();
        if (!pk) pk = malloc(crypto_kem_PUBLICKEYBYTES);
        if (!pk) abort();
        if (!sk) sk = malloc(crypto_kem_SECRETKEYBYTES);
        if (!sk) abort();

        randombytes_init(seed[i], NULL, 256);

        fprintf(fp_rsp, "count = %d\n", i);
        fprintBstr(fp_rsp, "seed = ", seed[i], 48);
       
        gettimeofday(&start_keygen, NULL);

        CALLGRIND_START_INSTRUMENTATION;
        ret_val = crypto_kem_keypair(pk, sk);
        CALLGRIND_STOP_INSTRUMENTATION;
        CALLGRIND_DUMP_STATS;

        gettimeofday(&end_keygen, 0);
        long seconds_keygen = end_keygen.tv_sec - start_keygen.tv_sec;
        long microseconds_keygen = end_keygen.tv_usec - start_keygen.tv_usec;
        double elapsed_keygen = seconds_keygen + microseconds_keygen*0.000001;
        sum_keygen += elapsed_keygen;
        times_keygen = times_keygen + 1;

        if (ret_val != 0) {
        // if ( (ret_val = crypto_kem_keypair(pk, sk)) != 0) {
            fprintf(stderr, "crypto_kem_keypair returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "pk = ", pk, crypto_kem_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk = ", sk, crypto_kem_SECRETKEYBYTES);
        

        // for(int t=0;t<100; t++){

            gettimeofday(&start_enc, NULL);

            //CALLGRIND_START_INSTRUMENTATION;
            ret_val = crypto_kem_enc(ct, ss, pk);
            //CALLGRIND_STOP_INSTRUMENTATION;
	    //CALLGRIND_DUMP_STATS;

            gettimeofday(&end_enc, 0);
            long seconds_enc = end_enc.tv_sec - start_enc.tv_sec;
            long microseconds_enc = end_enc.tv_usec - start_enc.tv_usec;
            double elapsed_enc = seconds_enc + microseconds_enc*0.000001;
            sum_enc += elapsed_enc;
            times_enc = times_enc + 1;

            if (ret_val != 0) {
            // if ( (ret_val = crypto_kem_enc(ct, ss, pk)) != 0) {
                fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
                return KAT_CRYPTO_FAILURE;
            }
            fprintBstr(fp_rsp, "ct = ", ct, crypto_kem_CIPHERTEXTBYTES);
            fprintBstr(fp_rsp, "ss = ", ss, crypto_kem_BYTES);
            
            fprintf(fp_rsp, "\n");
    
            gettimeofday(&start_dec, NULL);

	    //CALLGRIND_START_INSTRUMENTATION;
            ret_val =  crypto_kem_dec(ss1, ct, sk);
            //CALLGRIND_STOP_INSTRUMENTATION;
            //CALLGRIND_DUMP_STATS;



            gettimeofday(&end_dec, 0);
            long seconds_dec = end_dec.tv_sec - start_dec.tv_sec;
            long microseconds_dec = end_dec.tv_usec - start_dec.tv_usec;
            double elapsed_dec = seconds_dec + microseconds_dec*0.000001;
            sum_dec += elapsed_dec;
            times_dec = times_dec + 1;

            if (ret_val != 0) {
    //        if ( (ret_val = crypto_kem_dec(ss1, ct, sk)) != 0) {
                fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
                return KAT_CRYPTO_FAILURE;
            }
            
            if ( memcmp(ss, ss1, crypto_kem_BYTES) ) {
                fprintf(stderr, "crypto_kem_dec returned bad 'ss' value\n");
                return KAT_CRYPTO_FAILURE;
            }

        // }//t
    }
	
    printf("\n\t**********TIMING RESULTS**********\t\n");    
    printf("Elim kernel :Avg Execution time is: %0.3f miliseconds \n",(sum_elim)*1000/times_elim);
    printf("Synd kernel :Avg Execution time is: %0.3f miliseconds \n",(sum_synd)*1000/times_synd);
    printf("Syndrome kernel :Avg Execution time is: %0.3f miliseconds \n",(sum_syndrome)*1000/(times_syndrome));
    
    printf("\nKeygen :Avg Execution time is: %0.3f miliseconds \n",(sum_keygen)*1000/times_keygen);
    printf("Enc :Avg Execution time is: %0.3f miliseconds \n",(sum_enc)*1000/times_enc);
    printf("Dec :Avg Execution time is: %0.3f miliseconds \n",(sum_dec)*1000/times_dec);
    printf("Encrypt :Avg Execution time is: %0.3f miliseconds \n",(sum_encrypt)*1000/times_encrypt);
    printf("Decrypt :Avg Execution time is: %0.3f miliseconds \n",(sum_decrypt)*1000/times_decrypt);

    return KAT_SUCCESS;
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long i;

	fprintf(fp, "%s", S);

	for ( i=0; i<L; i++ )
		fprintf(fp, "%02X", A[i]);

	if ( L == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}
