#include "vec3.h"
#include "utils.h"
#include <math.h>

Vec3 *vec3_new(float x0, float x1, float x2) new(Vec3, x0, x1, x2);
Vec3 vec3_neg(Vec3 u) { return (Vec3){-u.x[0], -u.x[1], -u.x[2]}; }
Vec3 vec3vec3_add(Vec3 u, Vec3 v) { return (Vec3){u.x[0] + v.x[0], u.x[1] + v.x[1], u.x[2] + v.x[2]}; }
Vec3 vec3vec3_sub(Vec3 u, Vec3 v) { return (Vec3){u.x[0] - v.x[0], u.x[1] - v.x[1], u.x[2] - v.x[2]}; }
Vec3 vec3vec3_mul(Vec3 u, Vec3 v) { return (Vec3){u.x[0] * v.x[0], u.x[1] * v.x[1], u.x[2] * v.x[2]}; }
Vec3 vec3vec3_div(Vec3 u, Vec3 v) { return (Vec3){u.x[0] / v.x[0], u.x[1] / v.x[1], u.x[2] / v.x[2]}; }
Vec3 vec3float_add(Vec3 u, float v) { return (Vec3){u.x[0] + v, u.x[1] + v, u.x[2] + v}; }
Vec3 vec3float_sub(Vec3 u, float v) { return (Vec3){u.x[0] - v, u.x[1] - v, u.x[2] - v}; }
Vec3 vec3float_mul(Vec3 u, float v) { return (Vec3){u.x[0] * v, u.x[1] * v, u.x[2] * v}; }
Vec3 vec3float_div(Vec3 u, float v) { return vec3_mul(u, 1.0f / v); }
Vec3 vec3_lerp(Vec3 a, Vec3 b, float w) { return vec3_add(vec3_mul(a, 1.0f - w), vec3_mul(b, w)); }

float vec3_length2(Vec3 u) { return vec3_dot(u, u); }
float vec3_length(Vec3 u) { return sqrtf(vec3_length2(u)); }
float vec3_dot(Vec3 u, Vec3 v) { return u.x[0] * v.x[0] + u.x[1] * v.x[1] + u.x[2] * v.x[2]; }
Vec3 vec3_cross(Vec3 u, Vec3 v) {
  return (Vec3){
      u.x[1] * v.x[2] - v.x[1] * u.x[2],
      u.x[2] * v.x[0] - v.x[2] * u.x[0],
      u.x[0] * v.x[1] - v.x[0] * u.x[1],
  };
}
Vec3 vec3_unit(Vec3 u) { return vec3_div(u, vec3_length(u)); }
bool vec3_near_zero(Vec3 u) { return (fabsf(u.x[0]) < 1e-8f) && (fabsf(u.x[1]) < 1e-8f) && (fabsf(u.x[2]) < 1e-8f); }

Vec3 vec3_rand(PCG32State *rng) { return (Vec3){pcg32_f32(rng), pcg32_f32(rng), pcg32_f32(rng)}; }
Vec3 vec3_rand_between(PCG32State *rng, float lo, float hi) {
  return (Vec3){pcg32_f32_between(rng, lo, hi), pcg32_f32_between(rng, lo, hi), pcg32_f32_between(rng, lo, hi)};
}
Vec3 vec3_rand_unit_vector(PCG32State *rng) {
  for (;;) {
    Vec3 v = vec3_rand_between(rng, -1.0f, 1.0f);
    float length2 = vec3_length2(v);
    if (length2 < 1)
      return vec3_div(v, sqrtf(length2));
  }
}
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32State *rng) {
  Vec3 v = vec3_rand_unit_vector(rng);
  return (vec3_dot(v, normal) > 0.0f) ? v : vec3_neg(v);
}
