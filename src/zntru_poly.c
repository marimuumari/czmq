#include "czmq_classes.h"

/* From the original poly_mod.c */
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

void poly_mod_3_Phi_n(poly *r)
{
  int i;
  for(i=0; i <NTRU_N; i++)
    r->coeffs[i] = mod3(r->coeffs[i] + 2*r->coeffs[NTRU_N-1]);
}

void poly_mod_q_Phi_n(poly *r)
{
  int i;
  for(i=0; i<NTRU_N; i++)
    r->coeffs[i] = r->coeffs[i] - r->coeffs[NTRU_N-1];
}

void poly_Rq_to_S3(poly *r, const poly *a)
{
  int i;
  uint16_t flag;

  /* The coefficients of a are stored as non-negative integers. */
  /* We must translate to representatives in [-q/2, q/2) before */
  /* reduction mod 3.                                           */
  for(i=0; i<NTRU_N; i++)
  {
    /* Need an explicit reduction mod q here                    */
    r->coeffs[i] = MODQ(a->coeffs[i]);

    /* flag = 1 if r[i] >= q/2 else 0                            */
    flag = r->coeffs[i] >> (NTRU_LOGQ-1);

    /* Now we will add (-q) mod 3 if r[i] >= q/2                 */
    /* Note (-q) mod 3=(-2^k) mod 3=1<<(1-(k&1))                */
    r->coeffs[i] += flag << (1-(NTRU_LOGQ&1));
  }

  poly_mod_3_Phi_n(r);
}



/* From the original poly_r2_inv.c*/
/* return -1 if x<0 and y<0; otherwise return 0 */
static inline int16_t both_negative_mask(int16_t x,int16_t y)
{
  return (x & y) >> 15;
}

void poly_R2_inv(poly *r, const poly *a)
{
  poly f, g, v, w;
  size_t i, loop;
  int16_t delta,sign,swap,t;

  for (i = 0;i < NTRU_N;++i) v.coeffs[i] = 0;
  for (i = 0;i < NTRU_N;++i) w.coeffs[i] = 0;
  w.coeffs[0] = 1;

  for (i = 0;i < NTRU_N;++i) f.coeffs[i] = 1;
  for (i = 0;i < NTRU_N-1;++i) g.coeffs[NTRU_N-2-i] = (a->coeffs[i] ^ a->coeffs[NTRU_N-1]) & 1;
  g.coeffs[NTRU_N-1] = 0;

  delta = 1;

  for (loop = 0;loop < 2*(NTRU_N-1)-1;++loop) {
    for (i = NTRU_N-1;i > 0;--i) v.coeffs[i] = v.coeffs[i-1];
    v.coeffs[0] = 0;

    sign = g.coeffs[0] & f.coeffs[0];
    swap = both_negative_mask(-delta,-(int16_t) g.coeffs[0]);
    delta ^= swap & (delta ^ -delta);
    delta += 1;

    for (i = 0;i < NTRU_N;++i) {
      t = swap&(f.coeffs[i]^g.coeffs[i]); f.coeffs[i] ^= t; g.coeffs[i] ^= t;
      t = swap&(v.coeffs[i]^w.coeffs[i]); v.coeffs[i] ^= t; w.coeffs[i] ^= t;
    }

    for (i = 0;i < NTRU_N;++i) g.coeffs[i] = g.coeffs[i]^(sign&f.coeffs[i]);
    for (i = 0;i < NTRU_N;++i) w.coeffs[i] = w.coeffs[i]^(sign&v.coeffs[i]);
    for (i = 0;i < NTRU_N-1;++i) g.coeffs[i] = g.coeffs[i+1];
    g.coeffs[NTRU_N-1] = 0;
  }

  for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = v.coeffs[NTRU_N-2-i];
  r->coeffs[NTRU_N-1] = 0;
}



/* From the original poly_rq_mul.c */
void poly_Rq_mul(poly *r, const poly *a, const poly *b)
{
  int k,i;

  for(k=0; k<NTRU_N; k++)
  {
    r->coeffs[k] = 0;
    for(i=1; i<NTRU_N-k; i++)
      r->coeffs[k] += a->coeffs[k+i] * b->coeffs[NTRU_N-i];
    for(i=0; i<k+1; i++)
      r->coeffs[k] += a->coeffs[k-i] * b->coeffs[i];
  }
}



/* From the original poly_s3_inv.c*/
static inline uint8_t mod4(uint8_t a) /* a between 0 and 9 */
{
  int16_t t, c;
  a = (a >> 2) + (a & 3); /* between 0 and 4 */
  t = a - 3;
  c = t >> 5;
  return (uint8_t) (t^(c&(a^t)));
}


void poly_S3_inv(poly *r, const poly *a)
{
  poly f, g, v, w;
  size_t i,loop;
  int16_t delta,sign,swap,t;

  for (i = 0;i < NTRU_N;++i) v.coeffs[i] = 0;
  for (i = 0;i < NTRU_N;++i) w.coeffs[i] = 0;
  w.coeffs[0] = 1;

  for (i = 0;i < NTRU_N;++i) f.coeffs[i] = 1;
  for (i = 0;i < NTRU_N-1;++i) g.coeffs[NTRU_N-2-i] = mod4((a->coeffs[i] & 3) + 2*(a->coeffs[NTRU_N-1] & 3));
  g.coeffs[NTRU_N-1] = 0;

  delta = 1;

  for (loop = 0;loop < 2*(NTRU_N-1)-1;++loop) {
    for (i = NTRU_N-1;i > 0;--i) v.coeffs[i] = v.coeffs[i-1];
    v.coeffs[0] = 0;

    sign = mod4((uint8_t) (2 * g.coeffs[0] * f.coeffs[0]));
    swap = both_negative_mask(-delta,-(int16_t) g.coeffs[0]);
    delta ^= swap & (delta ^ -delta);
    delta += 1;

    for (i = 0;i < NTRU_N;++i) {
      t = swap&(f.coeffs[i]^g.coeffs[i]); f.coeffs[i] ^= t; g.coeffs[i] ^= t;
      t = swap&(v.coeffs[i]^w.coeffs[i]); v.coeffs[i] ^= t; w.coeffs[i] ^= t;
    }

    for (i = 0;i < NTRU_N;++i) g.coeffs[i] = mod4((uint8_t) (g.coeffs[i]+sign*f.coeffs[i]));
    for (i = 0;i < NTRU_N;++i) w.coeffs[i] = mod4((uint8_t) (w.coeffs[i]+sign*v.coeffs[i]));
    for (i = 0;i < NTRU_N-1;++i) g.coeffs[i] = g.coeffs[i+1];
    g.coeffs[NTRU_N-1] = 0;
  }

  sign = f.coeffs[0];
  for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = mod4((uint8_t) (sign*v.coeffs[NTRU_N-2-i]));
  r->coeffs[NTRU_N-1] = 0;
}



/*From the original poly.c */
/* Map {0, 1, 2} -> {0,1,q-1} in place */
void poly_Z3_to_Zq(poly *r)
{
  int i;
  for(i=0; i<NTRU_N; i++)
    r->coeffs[i] = r->coeffs[i] | ((-(r->coeffs[i]>>1)) & (NTRU_Q-1));
}

/* Map {0, 1, q-1} -> {0,1,2} in place */
void poly_trinary_Zq_to_Z3(poly *r)
{
  int i;
  for(i=0; i<NTRU_N; i++)
  {
    r->coeffs[i] = MODQ(r->coeffs[i]);
    r->coeffs[i] = 3 & (r->coeffs[i] ^ (r->coeffs[i]>>(NTRU_LOGQ-1)));
  }
}

void poly_Sq_mul(poly *r, const poly *a, const poly *b)
{
  poly_Rq_mul(r, a, b);
  poly_mod_q_Phi_n(r);
}

void poly_S3_mul(poly *r, const poly *a, const poly *b)
{
  poly_Rq_mul(r, a, b);
  poly_mod_3_Phi_n(r);
}

static void poly_R2_inv_to_Rq_inv(poly *r, const poly *ai, const poly *a)
{
#if NTRU_Q <= 256 || NTRU_Q >= 65536
#error "poly_R2_inv_to_Rq_inv in poly.c assumes 256 < q < 65536"
#endif

  int i;
  poly b, c;
  poly s;

  // for 0..4
  //    ai = ai * (2 - a*ai)  mod q
  for(i=0; i<NTRU_N; i++)
    b.coeffs[i] = -(a->coeffs[i]);

  for(i=0; i<NTRU_N; i++)
    r->coeffs[i] = ai->coeffs[i];

  poly_Rq_mul(&c, r, &b);
  c.coeffs[0] += 2; // c = 2 - a*ai
  poly_Rq_mul(&s, &c, r); // s = ai*c

  poly_Rq_mul(&c, &s, &b);
  c.coeffs[0] += 2; // c = 2 - a*s
  poly_Rq_mul(r, &c, &s); // r = s*c

  poly_Rq_mul(&c, r, &b);
  c.coeffs[0] += 2; // c = 2 - a*r
  poly_Rq_mul(&s, &c, r); // s = r*c

  poly_Rq_mul(&c, &s, &b);
  c.coeffs[0] += 2; // c = 2 - a*s
  poly_Rq_mul(r, &c, &s); // r = s*c
}

void poly_Rq_inv(poly *r, const poly *a)
{
  poly ai2;
  poly_R2_inv(&ai2, a);
  poly_R2_inv_to_Rq_inv(r, &ai2, a);
}



/*From the original pack3.c*/
void poly_S3_tobytes(unsigned char msg[NTRU_OWCPA_MSGBYTES], const poly *a)
{
  int i;
  unsigned char c;
#if NTRU_PACK_DEG > (NTRU_PACK_DEG / 5) * 5  // if 5 does not divide NTRU_N-1
  int j;
#endif

  for(i=0; i<NTRU_PACK_DEG/5; i++)
  {
    c =        a->coeffs[5*i+4] & 255;
    c = (3*c + a->coeffs[5*i+3]) & 255;
    c = (3*c + a->coeffs[5*i+2]) & 255;
    c = (3*c + a->coeffs[5*i+1]) & 255;
    c = (3*c + a->coeffs[5*i+0]) & 255;
    msg[i] = c;
  }
#if NTRU_PACK_DEG > (NTRU_PACK_DEG / 5) * 5  // if 5 does not divide NTRU_N-1
  i = NTRU_PACK_DEG/5;
  c = 0;
  for(j = NTRU_PACK_DEG - (5*i) - 1; j>=0; j--)
    c = (3*c + a->coeffs[5*i+j]) & 255;
  msg[i] = c;
#endif
}

void poly_S3_frombytes(poly *r, const unsigned char msg[NTRU_OWCPA_MSGBYTES])
{
  int i;
  unsigned char c;
#if NTRU_PACK_DEG > (NTRU_PACK_DEG / 5) * 5  // if 5 does not divide NTRU_N-1
  int j;
#endif

  for(i=0; i<NTRU_PACK_DEG/5; i++)
  {
    c = msg[i];
    r->coeffs[5*i+0] = c;
    r->coeffs[5*i+1] = c * 171 >> 9;  // this is division by 3
    r->coeffs[5*i+2] = c * 57 >> 9;  // division by 3^2
    r->coeffs[5*i+3] = c * 19 >> 9;  // division by 3^3
    r->coeffs[5*i+4] = c * 203 >> 14;  // etc.
  }
#if NTRU_PACK_DEG > (NTRU_PACK_DEG / 5) * 5  // if 5 does not divide NTRU_N-1
  i = NTRU_PACK_DEG/5;
  c = msg[i];
  for(j=0; (5*i+j)<NTRU_PACK_DEG; j++)
  {
    r->coeffs[5*i+j] = c;
    c = c * 171 >> 9;
  }
#endif
  r->coeffs[NTRU_N-1] = 0;
  poly_mod_3_Phi_n(r);
}



/* From the original packq.c */
void poly_Sq_tobytes(unsigned char *r, const poly *a)
{
  int i,j;
  uint16_t t[8];

  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    for(j=0;j<8;j++)
      t[j] = MODQ(a->coeffs[8*i+j]);

    r[11 * i + 0] = (unsigned char) ( t[0]        & 0xff);
    r[11 * i + 1] = (unsigned char) ((t[0] >>  8) | ((t[1] & 0x1f) << 3));
    r[11 * i + 2] = (unsigned char) ((t[1] >>  5) | ((t[2] & 0x03) << 6));
    r[11 * i + 3] = (unsigned char) ((t[2] >>  2) & 0xff);
    r[11 * i + 4] = (unsigned char) ((t[2] >> 10) | ((t[3] & 0x7f) << 1));
    r[11 * i + 5] = (unsigned char) ((t[3] >>  7) | ((t[4] & 0x0f) << 4));
    r[11 * i + 6] = (unsigned char) ((t[4] >>  4) | ((t[5] & 0x01) << 7));
    r[11 * i + 7] = (unsigned char) ((t[5] >>  1) & 0xff);
    r[11 * i + 8] = (unsigned char) ((t[5] >>  9) | ((t[6] & 0x3f) << 2));
    r[11 * i + 9] = (unsigned char) ((t[6] >>  6) | ((t[7] & 0x07) << 5));
    r[11 * i + 10] = (unsigned char) ((t[7] >>  3));
  }

  for(j=0;j<NTRU_PACK_DEG-8*i;j++)
    t[j] = MODQ(a->coeffs[8*i+j]);
  for(; j<8; j++)
    t[j] = 0;

  switch(NTRU_PACK_DEG&0x07)
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r[11 * i + 0] = (unsigned char) (t[0]        & 0xff);
      r[11 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x1f) << 3);
      r[11 * i + 2] = (unsigned char) (t[1] >>  5) | ((t[2] & 0x03) << 6);
      r[11 * i + 3] = (unsigned char) (t[2] >>  2) & 0xff;
      r[11 * i + 4] = (unsigned char) (t[2] >> 10) | ((t[3] & 0x7f) << 1);
      r[11 * i + 5] = (unsigned char) (t[3] >>  7) | ((t[4] & 0x0f) << 4);
      break;
    case 2:
      r[11 * i + 0] = (unsigned char) (t[0]        & 0xff);
      r[11 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x1f) << 3);
      r[11 * i + 2] = (unsigned char) (t[1] >>  5) | ((t[2] & 0x03) << 6);
      break;
  }
}

void poly_Sq_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
    r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
    r->coeffs[8*i+2] = (a[11*i+ 2] >> 6) | (((uint16_t)a[11*i+ 3] & 0xff) << 2) | (((uint16_t)a[11*i+ 4] & 0x01) << 10);
    r->coeffs[8*i+3] = (a[11*i+ 4] >> 1) | (((uint16_t)a[11*i+ 5] & 0x0f) << 7);
    r->coeffs[8*i+4] = (a[11*i+ 5] >> 4) | (((uint16_t)a[11*i+ 6] & 0x7f) << 4);
    r->coeffs[8*i+5] = (a[11*i+ 6] >> 7) | (((uint16_t)a[11*i+ 7] & 0xff) << 1) | (((uint16_t)a[11*i+ 8] & 0x03) <<  9);
    r->coeffs[8*i+6] = (a[11*i+ 8] >> 2) | (((uint16_t)a[11*i+ 9] & 0x1f) << 6);
    r->coeffs[8*i+7] = (a[11*i+ 9] >> 5) | (((uint16_t)a[11*i+10] & 0xff) << 3);
  }
  switch(NTRU_PACK_DEG&0x07)
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
      r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
      r->coeffs[8*i+2] = (a[11*i+ 2] >> 6) | (((uint16_t)a[11*i+ 3] & 0xff) << 2) | (((uint16_t)a[11*i+ 4] & 0x01) << 10);
      r->coeffs[8*i+3] = (a[11*i+ 4] >> 1) | (((uint16_t)a[11*i+ 5] & 0x0f) << 7);
      break;
    case 2:
      r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
      r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
      break;
  }
  r->coeffs[NTRU_N-1] = 0;
}

void poly_Rq_sum_zero_tobytes(unsigned char *r, const poly *a)
{
  poly_Sq_tobytes(r, a);
}

void poly_Rq_sum_zero_frombytes(poly *r, const unsigned char *a)
{
  int i;
  poly_Sq_frombytes(r,a);

  /* Set r[n-1] so that the sum of coefficients is zero mod q */
  r->coeffs[NTRU_N-1] = 0;
  for(i=0;i<NTRU_PACK_DEG;i++)
    r->coeffs[NTRU_N-1] -= r->coeffs[i];
}



/*From the original poly_lift.c */
#ifdef NTRU_HPS
void poly_lift(poly *r, const poly *a)
{
  int i;
  for(i=0; i<NTRU_N; i++) {
    r->coeffs[i] = a->coeffs[i];
  }
  poly_Z3_to_Zq(r);
}
#endif

#ifdef NTRU_HRSS
void poly_lift(poly *r, const poly *a)
{
  /* NOTE: Assumes input is in {0,1,2}^N */
  /*       Produces output in [0,Q-1]^N */
  int i;
  poly b;
  uint16_t t, zj;

  /* Define z by <z*x^i, x-1> = delta_{i,0} mod 3:      */
  /*   t      = -1/N mod p = -N mod 3                   */
  /*   z[0]   = 2 - t mod 3                             */
  /*   z[1]   = 0 mod 3                                 */
  /*   z[j]   = z[j-1] + t mod 3                        */
  /* We'll compute b = a/(x-1) mod (3, Phi) using       */
  /*   b[0] = <z, a>, b[1] = <z*x,a>, b[2] = <z*x^2,a>  */
  /*   b[i] = b[i-3] - (a[i] + a[i-1] + a[i-2])         */
  t = 3 - (NTRU_N % 3);
  b.coeffs[0] = a->coeffs[0] * (2-t) + a->coeffs[1] * 0 + a->coeffs[2] * t;
  b.coeffs[1] = a->coeffs[1] * (2-t) + a->coeffs[2] * 0;
  b.coeffs[2] = a->coeffs[2] * (2-t);

  zj = 0; /* z[1] */
  for(i=3; i<NTRU_N; i++)
  {
    b.coeffs[0] += a->coeffs[i] * (zj + 2*t);
    b.coeffs[1] += a->coeffs[i] * (zj + t);
    b.coeffs[2] += a->coeffs[i] * zj;
    zj = (zj + t) % 3;
  }
  b.coeffs[1] += a->coeffs[0] * (zj + t);
  b.coeffs[2] += a->coeffs[0] * zj;
  b.coeffs[2] += a->coeffs[1] * (zj + t);

  b.coeffs[0] = b.coeffs[0];
  b.coeffs[1] = b.coeffs[1];
  b.coeffs[2] = b.coeffs[2];
  for(i=3; i<NTRU_N; i++) {
    b.coeffs[i] = b.coeffs[i-3] + 2*(a->coeffs[i] + a->coeffs[i-1] + a->coeffs[i-2]);
  }

  /* Finish reduction mod Phi by subtracting Phi * b[N-1] */
  poly_mod_3_Phi_n(&b);

  /* Switch from {0,1,2} to {0,1,q-1} coefficient representation */
  poly_Z3_to_Zq(&b);

  /* Multiply by (x-1) */
  r->coeffs[0] = -(b.coeffs[0]);
  for(i=0; i<NTRU_N-1; i++) {
    r->coeffs[i+1] = b.coeffs[i] - b.coeffs[i+1];
  }
}
#endif

