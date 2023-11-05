// https://github.com/imneme/pcg-c-basic

#ifndef PCG32_H
#define PCG32_H
#include <stdint.h>

typedef struct PCG32State {
  uint64_t state;
  uint64_t inc;
} PCG32State;

void pcg32_srandom_r(PCG32State *rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_random_r(PCG32State *rng);
float pcg32_randomf_r(PCG32State *rng);

#endif // PCG32_H
