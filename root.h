/*
  This file is for evaluating a polynomial at one or more field elements
*/

#ifndef ROOT_H
#define ROOT_H
#define eval CRYPTO_NAMESPACE(eval)
#define root CRYPTO_NAMESPACE(root)

#include "gf.h"

//void eval_sw_host(gf *, gf *, gf *);
//void eval_host(gf *, gf *, gf *);
gf eval(gf *, gf);
void root(gf *, gf *, gf *);

extern double sum_eval;
extern int times_eval;

#endif

