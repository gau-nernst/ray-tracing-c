#include <stdbool.h>
#include <stdint.h>

typedef struct Vec3 {
  float x;
  float y;
  float z;
} Vec3;

// https://stackoverflow.com/a/11763277
#define _GET_MACRO(_1, _2, _3, _4, FUNC, ...) FUNC
#define vec3_add2(x, y) _Generic((y), Vec3: vec3vec3_add, float: vec3float_add)(x, y)
#define vec3_add3(x, y, z) vec3_add2(vec3_add2(x, y), z)
#define vec3_add4(x, y, z, t) vec3_add2(vec3_add3(x, y, z), t)
#define vec3_add(...) _GET_MACRO(__VA_ARGS__, vec3_add4, vec3_add3, vec3_add2)(__VA_ARGS__)

#define vec3_sub(x, y) _Generic((y), Vec3: vec3vec3_sub, float: vec3float_sub)(x, y)
#define vec3_mul(x, y) _Generic((y), Vec3: vec3vec3_mul, float: vec3float_mul)(x, y)

Vec3 *vec3vec3_add_(Vec3 *self, Vec3 other);
Vec3 *vec3vec3_sub_(Vec3 *self, Vec3 other);
Vec3 *vec3float_mul_(Vec3 *self, float other);

Vec3 vec3_neg(Vec3 u);
Vec3 vec3vec3_add(Vec3 u, Vec3 v);
Vec3 vec3vec3_sub(Vec3 u, Vec3 v);
Vec3 vec3vec3_mul(Vec3 u, Vec3 v);
Vec3 vec3float_add(Vec3 u, float v);
Vec3 vec3float_sub(Vec3 u, float v);
Vec3 vec3float_mul(Vec3 u, float v);

float vec3_length2(Vec3 u);
float vec3_length(Vec3 u);
float vec3_dot(Vec3 u, Vec3 v);
Vec3 vec3_cross(Vec3 u, Vec3 v);
Vec3 vec3_unit(Vec3 u);

typedef struct Ray {
  Vec3 origin;
  Vec3 direction;
} Ray;

Vec3 ray_at(Ray ray, float t);

typedef struct HitRecord {
  Vec3 p;
  Vec3 normal;
  float t;
  bool front_face;
} HitRecord;

typedef struct Sphere {
  Vec3 center;
  float radius;
} Sphere;

bool hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool hit_spheres(const Sphere *spheres, int n, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);

typedef struct PCG32State {
  uint64_t state;
  uint64_t inc;
} PCG32State;

void pcg32_srandom_r(PCG32State *rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_random_r(PCG32State *rng);
float pcg32_randomf_r(PCG32State *rng);
