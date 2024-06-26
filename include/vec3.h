#ifndef VEC3_H
#define VEC3_H

#include "pcg32.h"
#include <stdbool.h>

// https://stackoverflow.com/a/11763277
#define _GET_MACRO(_1, _2, _3, _4, FUNC, ...) FUNC
#define vec3_add2(x, y) _Generic((y), Vec3: vec3_add_vec3, float: vec3_add_float)(x, y)
#define vec3_add3(x, y, z) vec3_add2(vec3_add2(x, y), z)
#define vec3_add4(x, y, z, t) vec3_add2(vec3_add3(x, y, z), t)
#define vec3_add(...) _GET_MACRO(__VA_ARGS__, vec3_add4, vec3_add3, vec3_add2)(__VA_ARGS__)

#define vec3_mul2(x, y) _Generic((y), Vec3: vec3_mul_vec3, float: vec3_mul_float)(x, y)
#define vec3_mul3(x, y, z) vec3_mul2(vec3_mul2(x, y), z)
#define vec3_mul4(x, y, z, t) vec3_mul2(vec3_mul3(x, y, z), t)
#define vec3_mul(...) _GET_MACRO(__VA_ARGS__, vec3_mul4, vec3_mul3, vec3_mul2)(__VA_ARGS__)

#define vec3_sub(x, y) _Generic((y), Vec3: vec3_sub_vec3, float: vec3_sub_float)(x, y)
#define vec3_div(x, y) _Generic((y), Vec3: vec3_div_vec3, float: vec3_div_float)(x, y)

typedef union Vec3 {
  struct {
    float x;
    float y;
    float z;
  };
  float values[3];
} Vec3;

extern const Vec3 VEC3_ZERO;
Vec3 vec3(float x, float y, float z);
Vec3 *Vec3_new(float x, float y, float z);
Vec3 vec3_neg(Vec3 u);
Vec3 vec3_inv(Vec3 u);

Vec3 vec3_add_vec3(Vec3 u, Vec3 v);
Vec3 vec3_mul_vec3(Vec3 u, Vec3 v);
Vec3 vec3_sub_vec3(Vec3 u, Vec3 v);
Vec3 vec3_div_vec3(Vec3 u, Vec3 v);

Vec3 vec3_add_float(Vec3 u, float v);
Vec3 vec3_sub_float(Vec3 u, float v);
Vec3 vec3_mul_float(Vec3 u, float v);
Vec3 vec3_div_float(Vec3 u, float v);

Vec3 vec3_lerp(Vec3 a, Vec3 b, float w);
Vec3 vec3_min(Vec3 u, Vec3 v);
Vec3 vec3_max(Vec3 u, Vec3 v);

float vec3_length2(Vec3 u);
float vec3_length(Vec3 u);
float vec3_dot(Vec3 u, Vec3 v);
Vec3 vec3_cross(Vec3 u, Vec3 v);
Vec3 vec3_normalize(Vec3 u);
bool vec3_near_zero(Vec3 u);

Vec3 vec3_rand(PCG32 *rng);
Vec3 vec3_rand_between(PCG32 *rng, float lo, float hi);
Vec3 vec3_rand_unit_vector(PCG32 *rng);
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32 *rng);

#endif // VEC3_H
