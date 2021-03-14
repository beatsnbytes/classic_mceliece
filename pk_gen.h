/*
  This file is for public-key generation
*/

#ifndef PK_GEN_H
#define PK_GEN_H
#define pk_gen CRYPTO_NAMESPACE(pk_gen)

#include "gf.h"

int pk_gen(unsigned char *, unsigned char *, uint32_t *, int16_t *, unsigned char *, unsigned char *);

extern double sum_list_elim_tokern[1];
extern double sum_list_elim_tohost[1];
extern double sum_list_elim_kernel[1];
extern int times_elim;
extern int times_elim_tohost;
extern int times_elim_tokern;

extern double sum_pk_loop;
extern int times_pk_loop;


extern double sum_parallel;
extern int times_parallel;

#endif

