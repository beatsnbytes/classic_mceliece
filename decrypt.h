/*
  This file is for Nieddereiter decryption
*/

#ifndef DECRYPT_H
#define DECRYPT_H
#define decrypt CRYPTO_NAMESPACE(decrypt)

int decrypt(unsigned char *, const unsigned char *, const unsigned char *);

extern double sum_synd;
extern int times_synd;

#endif

