#include "vec3.h"
#include <math.h>

Vec3 vec3_neg(Vec3 u) { return (Vec3){-u.x, -u.y, -u.z}; }
Vec3 vec3vec3_add(Vec3 u, Vec3 v) { return (Vec3){u.x + v.x, u.y + v.y, u.z + v.z}; }
Vec3 vec3vec3_sub(Vec3 u, Vec3 v) { return (Vec3){u.x - v.x, u.y - v.y, u.z - v.z}; }
Vec3 vec3vec3_mul(Vec3 u, Vec3 v) { return (Vec3){u.x * v.x, u.y * v.y, u.z * v.z}; }
Vec3 vec3vec3_div(Vec3 u, Vec3 v) { return (Vec3){u.x / v.x, u.y / v.y, u.z / v.z}; }
Vec3 vec3float_add(Vec3 u, float v) { return (Vec3){u.x + v, u.y + v, u.z + v}; }
Vec3 vec3float_sub(Vec3 u, float v) { return (Vec3){u.x - v, u.y - v, u.z - v}; }
Vec3 vec3float_mul(Vec3 u, float v) { return (Vec3){u.x * v, u.y * v, u.z * v}; }
Vec3 vec3float_div(Vec3 u, float v) { return vec3_mul(u, 1.0f / v); }
Vec3 vec3_lerp(Vec3 a, Vec3 b, float w) { return vec3_add(vec3_mul(a, 1.0f - w), vec3_mul(b, w)); }

float vec3_length2(Vec3 u) { return vec3_dot(u, u); }
float vec3_length(Vec3 u) { return sqrtf(vec3_length2(u)); }
float vec3_dot(Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
Vec3 vec3_cross(Vec3 u, Vec3 v) { return (Vec3){u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y}; }
Vec3 vec3_unit(Vec3 u) { return vec3_div(u, vec3_length(u)); }
bool vec3_near_zero(Vec3 u) { return (fabsf(u.x) < 1e-8f) && (fabsf(u.y) < 1e-8f) && (fabsf(u.z) < 1e-8f); }

Vec3 vec3_rand(PCG32State *rng) { return (Vec3){pcg32_randomf_r(rng), pcg32_randomf_r(rng), pcg32_randomf_r(rng)}; }
Vec3 vec3_rand_between(float low, float high, PCG32State *rng) {
  return vec3_add(vec3_mul(vec3_rand(rng), high - low), low);
}
Vec3 vec3_rand_unit_vector(PCG32State *rng) {
  for (;;) {
    Vec3 v = vec3_rand_between(-1.0f, 1.0f, rng);
    float length2 = vec3_length2(v);
    if (length2 < 1)
      return vec3_div(v, sqrtf(length2));
  }
}
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32State *rng) {
  Vec3 v = vec3_rand_unit_vector(rng);
  return (vec3_dot(v, normal) > 0.0f) ? v : vec3_neg(v);
}
