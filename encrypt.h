/*
  This file is for Niederreiter encryption
*/

#ifndef ENCRYPT_H
#define ENCRYPT_H
#define encrypt CRYPTO_NAMESPACE(encrypt)

void encrypt(unsigned char *, const unsigned char *, unsigned char *);

extern double sum_syndrome;
extern int times_syndrome;

extern double sum_syndrome_2;
extern int times_syndrome_2;


#endif

