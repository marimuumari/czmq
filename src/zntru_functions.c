#include "czmq_classes.h"
#include "../include/zntru_poly.h"

#include <stdint.h>
#define int32 int32_t

/* From the original cmov.c */
/* b = 1 means mov, b = 0 means don't mov*/
void cmov(unsigned char *r, const unsigned char *x, size_t len, unsigned char b)
{
  size_t i;

  b = (~b + 1);
  for(i=0;i<len;i++)
    r[i] ^= b & (x[i] ^ r[i]);
}



/* From the original crypto_sort_int32.c */
#define int32_MINMAX(a,b) \
do { \
  int32_t ab = (b) ^ (a); \
  int32_t c = (int32_t)((int64_t)(b) - (int64_t)(a)); \
  c ^= ab & (c ^ (b)); \
  c >>= 31; \
  c &= ab; \
  (a) ^= c; \
  (b) ^= c; \
} while(0)

/* assume 2 <= n <= 0x40000000 */
void crypto_sort_int32(int32 *array,size_t n)
{
  size_t top,p,q,r,i,j;
  int32 *x = array;

  top = 1;
  while (top < n - top) top += top;

  for (p = top;p >= 1;p >>= 1) {
    i = 0;
    while (i + 2 * p <= n) {
      for (j = i;j < i + p;++j) {
        int32_MINMAX(x[j],x[j+p]);
      }
      i += 2 * p;
    }
    for (j = i;j < n - p;++j) {
      int32_MINMAX(x[j],x[j+p]);
    }

    i = 0;
    j = 0;
    for (q = top;q > p;q >>= 1) {
      if (j != i) {
        for (;;) {
          if (j == n - q) goto done;
          int32 a = x[j + p];
          for (r = q;r > p;r >>= 1) {
            int32_MINMAX(a,x[j + r]);
          }
          x[j + p] = a;
          ++j;
          if (j == i + p) {
            i += 2 * p;
            break;
          }
        }
      }
      while (i + p <= n - q) {
        for (j = i;j < i + p;++j) {
          int32 a = x[j + p];
          for (r = q;r > p;r >>= 1) {
            int32_MINMAX(a,x[j+r]);
          }
          x[j + p] = a;
        }
        i += 2 * p;
      }
      /* now i + p > n - q */
      j = i;
      while (j < n - q) {
        int32 a = x[j + p];
        for (r = q;r > p;r >>= 1) {
          int32_MINMAX(a,x[j+r]);
        }
        x[j + p] = a;
        ++j;
      }

      done: ;
    }
  }
}



/* From the original owcpa.c */
static int owcpa_check_ciphertext(const unsigned char *ciphertext)
{
  /* A ciphertext is log2(q)*(n-1) bits packed into bytes.  */
  /* Check that any unused bits of the final byte are zero. */

  uint16_t t = 0;

  t = ciphertext[NTRU_CIPHERTEXTBYTES-1];
  t &= 0xff << (8-(7 & (NTRU_LOGQ*NTRU_PACK_DEG)));

  /* We have 0 <= t < 256 */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 15));
}

static int owcpa_check_r(const poly *r)
{
  /* A valid r has coefficients in {0,1,q-1} and has r[N-1] = 0 */
  /* Note: We may assume that 0 <= r[i] <= q-1 for all i        */

  int i;
  uint32_t t = 0;
  uint16_t c;
  for(i=0; i<NTRU_N-1; i++)
  {
    c = r->coeffs[i];
    t |= (c + 1) & (NTRU_Q-4);  /* 0 iff c is in {-1,0,1,2} */
    t |= (c + 2) & 4;  /* 1 if c = 2, 0 if c is in {-1,0,1} */
  }
  t |= r->coeffs[NTRU_N-1]; /* Coefficient n-1 must be zero */

  /* We have 0 <= t < 2^16. */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 31));
}

#ifdef NTRU_HPS
static int owcpa_check_m(const poly *m)
{
  /* Check that m is in message space, i.e.                  */
  /*  (1)  |{i : m[i] = 1}| = |{i : m[i] = 2}|, and          */
  /*  (2)  |{i : m[i] != 0}| = NTRU_WEIGHT.                  */
  /* Note: We may assume that m has coefficients in {0,1,2}. */

  int i;
  uint32_t t = 0;
  uint16_t ps = 0;
  uint16_t ms = 0;
  for(i=0; i<NTRU_N; i++)
  {
    ps += m->coeffs[i] & 1;
    ms += m->coeffs[i] & 2;
  }
  t |= ps ^ (ms >> 1);   /* 0 if (1) holds */
  t |= ms ^ NTRU_WEIGHT; /* 0 if (1) and (2) hold */

  /* We have 0 <= t < 2^16. */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 31));
}
#endif

void owcpa_keypair(unsigned char *pk,
                   unsigned char *sk,
                   const unsigned char seed[NTRU_SAMPLE_FG_BYTES])
{
  int i;

  poly x1, x2, x3, x4, x5;

  poly *f=&x1, *g=&x2, *invf_mod3=&x3;
  poly *gf=&x3, *invgf=&x4, *tmp=&x5;
  poly *invh=&x3, *h=&x3;

  sample_fg(f,g,seed);

  poly_S3_inv(invf_mod3, f);
  poly_S3_tobytes(sk, f);
  poly_S3_tobytes(sk+NTRU_PACK_TRINARY_BYTES, invf_mod3);

  /* Lift coeffs of f and g from Z_p to Z_q */
  poly_Z3_to_Zq(f);
  poly_Z3_to_Zq(g);

#ifdef NTRU_HRSS
  /* g = 3*(x-1)*g */
  for(i=NTRU_N-1; i>0; i--)
    g->coeffs[i] = 3*(g->coeffs[i-1] - g->coeffs[i]);
  g->coeffs[0] = -(3*g->coeffs[0]);
#endif

#ifdef NTRU_HPS
  /* g = 3*g */
  for(i=0; i<NTRU_N; i++)
    g->coeffs[i] = 3 * g->coeffs[i];
#endif

  poly_Rq_mul(gf, g, f);

  poly_Rq_inv(invgf, gf);

  poly_Rq_mul(tmp, invgf, f);
  poly_Sq_mul(invh, tmp, f);
  poly_Sq_tobytes(sk+2*NTRU_PACK_TRINARY_BYTES, invh);

  poly_Rq_mul(tmp, invgf, g);
  poly_Rq_mul(h, tmp, g);
  poly_Rq_sum_zero_tobytes(pk, h);
}


void owcpa_enc(unsigned char *c,
               const poly *r,
               const poly *m,
               const unsigned char *pk)
{
  int i;
  poly x1, x2;
  poly *h = &x1, *liftm = &x1;
  poly *ct = &x2;

  poly_Rq_sum_zero_frombytes(h, pk);

  poly_Rq_mul(ct, r, h);

  poly_lift(liftm, m);
  for(i=0; i<NTRU_N; i++)
    ct->coeffs[i] = ct->coeffs[i] + liftm->coeffs[i];

  poly_Rq_sum_zero_tobytes(c, ct);
}

int owcpa_dec(unsigned char *rm,
              const unsigned char *ciphertext,
              const unsigned char *secretkey)
{
  int i;
  int fail;
  poly x1, x2, x3, x4;

  poly *c = &x1, *f = &x2, *cf = &x3;
  poly *mf = &x2, *finv3 = &x3, *m = &x4;
  poly *liftm = &x2, *invh = &x3, *r = &x4;
  poly *b = &x1;

  poly_Rq_sum_zero_frombytes(c, ciphertext);
  poly_S3_frombytes(f, secretkey);
  poly_Z3_to_Zq(f);

  poly_Rq_mul(cf, c, f);
  poly_Rq_to_S3(mf, cf);

  poly_S3_frombytes(finv3, secretkey+NTRU_PACK_TRINARY_BYTES);
  poly_S3_mul(m, mf, finv3);
  poly_S3_tobytes(rm+NTRU_PACK_TRINARY_BYTES, m);

  fail = 0;

  /* Check that the unused bits of the last byte of the ciphertext are zero */
  fail |= owcpa_check_ciphertext(ciphertext);

  /* For the IND-CCA2 KEM we must ensure that c = Enc(h, (r,m)).             */
  /* We can avoid re-computing r*h + Lift(m) as long as we check that        */
  /* r (defined as b/h mod (q, Phi_n)) and m are in the message space.       */
  /* (m can take any value in S3 in NTRU_HRSS) */
#ifdef NTRU_HPS
  fail |= owcpa_check_m(m);
#endif

  /* b = c - Lift(m) mod (q, x^n - 1) */
  poly_lift(liftm, m);
  for(i=0; i<NTRU_N; i++)
    b->coeffs[i] = c->coeffs[i] - liftm->coeffs[i];

  /* r = b / h mod (q, Phi_n) */
  poly_Sq_frombytes(invh, secretkey+2*NTRU_PACK_TRINARY_BYTES);
  poly_Sq_mul(r, b, invh);

  /* NOTE: Our definition of r as b/h mod (q, Phi_n) follows Figure 4 of     */
  /*   [Sch18] https://eprint.iacr.org/2018/1174/20181203:032458.            */
  /* This differs from Figure 10 of Saito--Xagawa--Yamakawa                  */
  /*   [SXY17] https://eprint.iacr.org/2017/1005/20180516:055500             */
  /* where r gets a final reduction modulo p.                                */
  /* We need this change to use Proposition 1 of [Sch18].                    */

  /* Proposition 1 of [Sch18] shows that re-encryption with (r,m) yields c.  */
  /* if and only if fail==0 after the following call to owcpa_check_r        */
  /* The procedure given in Fig. 8 of [Sch18] can be skipped because we have */
  /* c(1) = 0 due to the use of poly_Rq_sum_zero_{to,from}bytes.             */
  fail |= owcpa_check_r(r);

  poly_trinary_Zq_to_Z3(r);
  poly_S3_tobytes(rm, r);

  return fail;
}



/* From the original sample_iid.c */
static uint16_t mod3(uint16_t a)
{
  uint16_t r;
  int16_t t, c;

  r = (a >> 8) + (a & 0xff); // r mod 255 == a mod 255
  r = (r >> 4) + (r & 0xf); // r' mod 15 == r mod 15
  r = (r >> 2) + (r & 0x3); // r' mod 3 == r mod 3
  r = (r >> 2) + (r & 0x3); // r' mod 3 == r mod 3

  t = r - 3;
  c = t >> 15;

  return (c&r) ^ (~c&t);
}

void sample_iid(poly *r, const unsigned char uniformbytes[NTRU_SAMPLE_IID_BYTES])
{
  int i;
  /* {0,1,...,255} -> {0,1,2}; Pr[0] = 86/256, Pr[1] = Pr[-1] = 85/256 */
  for(i=0; i<NTRU_N-1; i++)
    r->coeffs[i] = mod3(uniformbytes[i]);

  r->coeffs[NTRU_N-1] = 0;
}



/* From the original sample.c */
void sample_fg(poly *f, poly *g, const unsigned char uniformbytes[NTRU_SAMPLE_FG_BYTES])
{
#ifdef NTRU_HRSS
  sample_iid_plus(f,uniformbytes);
  sample_iid_plus(g,uniformbytes+NTRU_SAMPLE_IID_BYTES);
#endif

#ifdef NTRU_HPS
  sample_iid(f,uniformbytes);
  sample_fixed_type(g,uniformbytes+NTRU_SAMPLE_IID_BYTES);
#endif
}

void sample_rm(poly *r, poly *m, const unsigned char uniformbytes[NTRU_SAMPLE_RM_BYTES])
{
#ifdef NTRU_HRSS
  sample_iid(r,uniformbytes);
  sample_iid(m,uniformbytes+NTRU_SAMPLE_IID_BYTES);
#endif

#ifdef NTRU_HPS
  sample_iid(r,uniformbytes);
  sample_fixed_type(m,uniformbytes+NTRU_SAMPLE_IID_BYTES);
#endif
}

#ifdef NTRU_HRSS
void sample_iid_plus(poly *r, const unsigned char uniformbytes[NTRU_SAMPLE_IID_BYTES])
{
  /* Sample r using sample_iid then conditionally flip    */
  /* signs of even index coefficients so that <x*r, r> >= 0.      */

  int i;
  uint16_t s = 0;

  sample_iid(r, uniformbytes);

  /* Map {0,1,2} -> {0, 1, 2^16 - 1} */
  for(i=0; i<NTRU_N-1; i++)
    r->coeffs[i] = r->coeffs[i] | (-(r->coeffs[i]>>1));

  /* s = <x*r, r>.  (r[n-1] = 0) */
  for(i=0; i<NTRU_N-1; i++)
    s += (uint16_t)((uint32_t)r->coeffs[i + 1] * (uint32_t)r->coeffs[i]);

  /* Extract sign of s (sign(0) = 1) */
  s = 1 | (-(s>>15));

  for(i=0; i<NTRU_N; i+=2)
    r->coeffs[i] = (uint16_t)((uint32_t)s * (uint32_t)r->coeffs[i]);

  /* Map {0,1,2^16-1} -> {0, 1, 2} */
  for(i=0; i<NTRU_N; i++)
    r->coeffs[i] = 3 & (r->coeffs[i] ^ (r->coeffs[i]>>15));
}
#endif

#ifdef NTRU_HPS
void sample_fixed_type(poly *r, const unsigned char u[NTRU_SAMPLE_FT_BYTES])
{
  // Assumes NTRU_SAMPLE_FT_BYTES = ceil(30*(n-1)/8)

  int32_t s[NTRU_N-1];
  int i;

  // Use 30 bits of u per word
  for (i = 0; i < (NTRU_N-1)/4; i++)
  {
    s[4*i+0] =                              (u[15*i+ 0] << 2) + (u[15*i+ 1] << 10) + (u[15*i+ 2] << 18) + ((uint32_t) u[15*i+ 3] << 26);
    s[4*i+1] = ((u[15*i+ 3] & 0xc0) >> 4) + (u[15*i+ 4] << 4) + (u[15*i+ 5] << 12) + (u[15*i+ 6] << 20) + ((uint32_t) u[15*i+ 7] << 28);
    s[4*i+2] = ((u[15*i+ 7] & 0xf0) >> 2) + (u[15*i+ 8] << 6) + (u[15*i+ 9] << 14) + (u[15*i+10] << 22) + ((uint32_t) u[15*i+11] << 30);
    s[4*i+3] =  (u[15*i+11] & 0xfc)       + (u[15*i+12] << 8) + (u[15*i+13] << 16) + ((uint32_t) u[15*i+14] << 24);
  }
#if (NTRU_N - 1) > ((NTRU_N - 1) / 4) * 4 // (N-1) = 2 mod 4
  i = (NTRU_N-1)/4;
  s[4*i+0] =                              (u[15*i+ 0] << 2) + (u[15*i+ 1] << 10) + (u[15*i+ 2] << 18) + ((uint32_t) u[15*i+ 3] << 26);
  s[4*i+1] = ((u[15*i+ 3] & 0xc0) >> 4) + (u[15*i+ 4] << 4) + (u[15*i+ 5] << 12) + (u[15*i+ 6] << 20) + ((uint32_t) u[15*i+ 7] << 28);
#endif

  for (i = 0; i<NTRU_WEIGHT/2; i++) s[i] |=  1;

  for (i = NTRU_WEIGHT/2; i<NTRU_WEIGHT; i++) s[i] |=  2;

  crypto_sort_int32(s,NTRU_N-1);

  for(i=0; i<NTRU_N-1; i++)
    r->coeffs[i] = ((uint16_t) (s[i] & 3));

  r->coeffs[NTRU_N-1] = 0;
}
#endif