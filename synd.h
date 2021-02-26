/*
  This file is for syndrome computation
*/

#ifndef SYND_H
#define SYND_H
#define synd CRYPTO_NAMESPACE(synd)

#include "gf.h"

void synd(gf *, gf *, gf *, unsigned char *);

extern double sum_list_synd_tokern[1];
extern double sum_list_synd_tohost[1];
extern double sum_list_synd_kernel[2];
extern int times_synd;
extern int times_synd_tohost;
extern int times_synd_tokern;

extern double sum_list_synd_last_tokern[1];
extern double sum_list_synd_last_tohost[1];
extern double sum_list_synd_last_kernel[1];
extern int times_synd_last;
extern int times_synd_last_tohost;
extern int times_synd_last_tokern;

#endif

