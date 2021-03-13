#ifndef crypto_kem_mceliece460896_H
#define crypto_kem_mceliece460896_H

#define crypto_kem_mceliece460896_vec_PUBLICKEYBYTES 524160
#define crypto_kem_mceliece460896_vec_SECRETKEYBYTES 13608
#define crypto_kem_mceliece460896_vec_CIPHERTEXTBYTES 188
#define crypto_kem_mceliece460896_vec_BYTES 32

 
#ifdef __cplusplus
extern "C" {
#endif
extern int crypto_kem_mceliece460896_vec_keypair(unsigned char *,unsigned char *);
extern int crypto_kem_mceliece460896_vec_enc(unsigned char *,unsigned char *,const unsigned char *);
extern int crypto_kem_mceliece460896_vec_dec(unsigned char *,const unsigned char *,const unsigned char *);
#ifdef __cplusplus
}
#endif

#define crypto_kem_mceliece460896_keypair crypto_kem_mceliece460896_vec_keypair
#define crypto_kem_mceliece460896_enc crypto_kem_mceliece460896_vec_enc
#define crypto_kem_mceliece460896_dec crypto_kem_mceliece460896_vec_dec
#define crypto_kem_mceliece460896_PUBLICKEYBYTES crypto_kem_mceliece460896_vec_PUBLICKEYBYTES
#define crypto_kem_mceliece460896_SECRETKEYBYTES crypto_kem_mceliece460896_vec_SECRETKEYBYTES
#define crypto_kem_mceliece460896_BYTES crypto_kem_mceliece460896_vec_BYTES
#define crypto_kem_mceliece460896_CIPHERTEXTBYTES crypto_kem_mceliece460896_vec_CIPHERTEXTBYTES

#endif
