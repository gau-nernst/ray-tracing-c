#include "pcg32.h"

void pcg32_seed(PCG32State *rng, uint64_t initstate, uint64_t initseq) {
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg32_u32(rng);
  rng->state += initstate;
  pcg32_u32(rng);
}

uint32_t pcg32_u32(PCG32State *rng) {
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + rng->inc;
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
uint32_t pcg32_u32_between(PCG32State *rng, uint32_t lo, uint32_t hi) {
  return lo + pcg32_u32(rng) % (hi - lo); // not accurate
}
float pcg32_f32(PCG32State *rng) { return (float)(pcg32_u32(rng) >> 8) / (float)(1 << 24); }
float pcg32_f32_between(PCG32State *rng, float lo, float hi) { return lo + pcg32_f32(rng) * (hi - lo); }
