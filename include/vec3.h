#ifndef VEC3_H
#define VEC3_H

#include "pcg32.h"
#include <stdbool.h>

// https://stackoverflow.com/a/11763277
#define _GET_MACRO(_1, _2, _3, _4, FUNC, ...) FUNC
#define vec3_add2(x, y) _Generic((y), Vec3: vec3vec3_add, float: vec3float_add)(x, y)
#define vec3_add3(x, y, z) vec3_add2(vec3_add2(x, y), z)
#define vec3_add4(x, y, z, t) vec3_add2(vec3_add3(x, y, z), t)
#define vec3_add(...) _GET_MACRO(__VA_ARGS__, vec3_add4, vec3_add3, vec3_add2)(__VA_ARGS__)

#define vec3_sub(x, y) _Generic((y), Vec3: vec3vec3_sub, float: vec3float_sub)(x, y)
#define vec3_mul(x, y) _Generic((y), Vec3: vec3vec3_mul, float: vec3float_mul)(x, y)
#define vec3_div(x, y) _Generic((y), Vec3: vec3vec3_div, float: vec3float_div)(x, y)

typedef struct Vec3 {
  float x;
  float y;
  float z;
} Vec3;

Vec3 vec3_neg(Vec3 u);
Vec3 vec3vec3_add(Vec3 u, Vec3 v);
Vec3 vec3vec3_sub(Vec3 u, Vec3 v);
Vec3 vec3vec3_mul(Vec3 u, Vec3 v);
Vec3 vec3vec3_div(Vec3 u, Vec3 v);
Vec3 vec3float_add(Vec3 u, float v);
Vec3 vec3float_sub(Vec3 u, float v);
Vec3 vec3float_mul(Vec3 u, float v);
Vec3 vec3float_div(Vec3 u, float v);
Vec3 vec3_lerp(Vec3 a, Vec3 b, float w);

float vec3_length2(Vec3 u);
float vec3_length(Vec3 u);
float vec3_dot(Vec3 u, Vec3 v);
Vec3 vec3_cross(Vec3 u, Vec3 v);
Vec3 vec3_unit(Vec3 u);
bool vec3_near_zero(Vec3 u);

Vec3 vec3_rand(PCG32State *rng);
Vec3 vec3_rand_between(PCG32State *rng, float lo, float hi);
Vec3 vec3_rand_unit_vector(PCG32State *rng);
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32State *rng);

#endif // VEC3_H
