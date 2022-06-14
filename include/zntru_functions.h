#ifndef ZNTRU_FUNCTIONS_H
#define ZNTRU_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif


#include "zntru_params.h"
#include "zntru_poly.h"

#include <stddef.h>
#include <stdint.h>


/* From the original cmov.h */
#define cmov CRYPTO_NAMESPACE(cmov)
void cmov(unsigned char *r, const unsigned char *x, size_t len, unsigned char b);


/* From the original crypto_sort_int32*/
#define crypto_sort_int32 CRYPTO_NAMESPACE(crypto_sort_int32)
void crypto_sort_int32(int32_t *array,size_t n);


/* From the original owcpa.h */
#define owcpa_samplemsg CRYPTO_NAMESPACE(owcpa_samplemsg)
void owcpa_samplemsg(unsigned char msg[NTRU_OWCPA_MSGBYTES],
                     const unsigned char seed[NTRU_SEEDBYTES]);

#define owcpa_keypair CRYPTO_NAMESPACE(owcpa_keypair)
void owcpa_keypair(unsigned char *pk,
                   unsigned char *sk,
                   const unsigned char seed[NTRU_SEEDBYTES]);

#define owcpa_enc CRYPTO_NAMESPACE(owcpa_enc)
void owcpa_enc(unsigned char *c,
               const poly *r,
               const poly *m,
               const unsigned char *pk);

#define owcpa_dec CRYPTO_NAMESPACE(owcpa_dec)
int owcpa_dec(unsigned char *rm,
              const unsigned char *ciphertext,
              const unsigned char *secretkey);


/* From the original sample.h */
#define sample_fg CRYPTO_NAMESPACE(sample_fg)
#define sample_rm CRYPTO_NAMESPACE(sample_rm)
void sample_fg(poly *f, poly *g, const unsigned char uniformbytes[NTRU_SAMPLE_FG_BYTES]);
void sample_rm(poly *r, poly *m, const unsigned char uniformbytes[NTRU_SAMPLE_RM_BYTES]);

#define sample_iid CRYPTO_NAMESPACE(sample_iid)
void sample_iid(poly *r, const unsigned char uniformbytes[NTRU_SAMPLE_IID_BYTES]);

#ifdef NTRU_HPS /* hps needs sample_fixed_type */
#define sample_fixed_type CRYPTO_NAMESPACE(sample_fixed_type)
void sample_fixed_type(poly *r, const unsigned char uniformbytes[NTRU_SAMPLE_FT_BYTES]);
#endif

#ifdef NTRU_HRSS /* hrss needs sample_iid_plus */
#define sample_iid_plus CRYPTO_NAMESPACE(sample_iid_plus)
void sample_iid_plus(poly *r, const unsigned char uniformbytes[NTRU_SAMPLE_IID_BYTES]);
#endif

#ifdef __cplusplus
}
#endif

#endif