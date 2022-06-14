#ifndef __ZNTRU_PARAMS_H_INCLUDED__
#define __ZNTRU_PARAMS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#define CRYPTO_SECRETKEYBYTES 935
#define CRYPTO_PUBLICKEYBYTES 699
#define CRYPTO_CIPHERTEXTBYTES 699
#define CRYPTO_BYTES 32

#define NTRU_HPS
#define NTRU_N 509
#define NTRU_LOGQ 11

#define CRYPTO_ALGNAME "ntruhps2048509"

#ifndef CRYPTO_NAMESPACE
#define CRYPTO_NAMESPACE(s) s
#endif

#define crypto_kem_keypair CRYPTO_NAMESPACE(keypair)
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

#define crypto_kem_enc CRYPTO_NAMESPACE(enc)
int crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);

#define crypto_kem_dec CRYPTO_NAMESPACE(dec)
int crypto_kem_dec(unsigned char *k, const unsigned char *c, const unsigned char *sk);



/* Do not modify below this line */

#define PAD32(X) ((((X) + 31)/32)*32)

#define NTRU_Q (1 << NTRU_LOGQ)
#define NTRU_WEIGHT (NTRU_Q/8 - 2)

#define NTRU_SEEDBYTES       32
#define NTRU_PRFKEYBYTES     32
#define NTRU_SHAREDKEYBYTES  32

#define NTRU_SAMPLE_IID_BYTES  (NTRU_N-1)
#define NTRU_SAMPLE_FT_BYTES   ((30*(NTRU_N-1)+7)/8)
#define NTRU_SAMPLE_FG_BYTES   (NTRU_SAMPLE_IID_BYTES+NTRU_SAMPLE_FT_BYTES)
#define NTRU_SAMPLE_RM_BYTES   (NTRU_SAMPLE_IID_BYTES+NTRU_SAMPLE_FT_BYTES)

#define NTRU_PACK_DEG (NTRU_N-1)
#define NTRU_PACK_TRINARY_BYTES    ((NTRU_PACK_DEG+4)/5)

#define NTRU_OWCPA_MSGBYTES       (2*NTRU_PACK_TRINARY_BYTES)
#define NTRU_OWCPA_PUBLICKEYBYTES ((NTRU_LOGQ*NTRU_PACK_DEG+7)/8)
#define NTRU_OWCPA_SECRETKEYBYTES (2*NTRU_PACK_TRINARY_BYTES + NTRU_OWCPA_PUBLICKEYBYTES)
#define NTRU_OWCPA_BYTES          ((NTRU_LOGQ*NTRU_PACK_DEG+7)/8)

#define NTRU_PUBLICKEYBYTES  (NTRU_OWCPA_PUBLICKEYBYTES)
#define NTRU_SECRETKEYBYTES  (NTRU_OWCPA_SECRETKEYBYTES + NTRU_PRFKEYBYTES)
#define NTRU_CIPHERTEXTBYTES (NTRU_OWCPA_BYTES)

#ifdef __cplusplus
}
#endif

#endif

