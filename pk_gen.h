/*
  This file is for public-key generation
*/

#ifndef PK_GEN_H
#define PK_GEN_H
#define pk_gen CRYPTO_NAMESPACE(pk_gen)

#include "gf.h"

int pk_gen(unsigned char *, unsigned char *, uint32_t *, int16_t *, unsigned char *, unsigned char *);

extern double sum_elim;
extern int times_elim;

extern double sum_elim_mig;
extern int times_elim_mig;

extern double sum_par;
extern int times_par;

#endif

