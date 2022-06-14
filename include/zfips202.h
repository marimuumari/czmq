#ifndef FIPS202_H
#define FIPS202_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define crypto_hash_sha3256 sha3_256
#define crypto_hash_sha3512 sha3_512
#define crypto_hash_shake256 shake256

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE  72

void shake128_absorb(uint64_t *s, const unsigned char *input, unsigned int inputByteLen);
void shake128_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s);

void shake256(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen);
void sha3_256(unsigned char *output, const unsigned char *input,  unsigned long long inlen);
void sha3_512(unsigned char *output, const unsigned char *input,  unsigned long long inlen);

#ifdef __cplusplus
}
#endif

#endif
