#ifndef PCG32_H
#define PCG32_H

#include <stdint.h>

// https://github.com/imneme/pcg-c-basic
typedef struct PCG32 {
  uint64_t state;
  uint64_t inc;
} PCG32;

void pcg32_seed(PCG32 *rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_u32(PCG32 *rng);
uint32_t pcg32_u32_between(PCG32 *rng, uint32_t lo, uint32_t hi);
float pcg32_f32(PCG32 *rng);
float pcg32_f32_between(PCG32 *rng, float lo, float hi);

#endif // PCG32_H
