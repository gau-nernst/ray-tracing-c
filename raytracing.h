#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#define clamp(x, lo, hi) min(max(x, lo), hi)

typedef struct PCG32State PCG32State;
typedef struct Vec3 Vec3;
typedef struct Ray Ray;
typedef struct HitRecord HitRecord;
typedef struct Sphere Sphere;
typedef enum MaterialType MaterialType;
typedef struct Material Material;
typedef struct World World;
typedef struct Camera Camera;
typedef enum TextureType TextureType;
typedef struct Texture Texture;

// https://stackoverflow.com/a/11763277
#define _GET_MACRO(_1, _2, _3, _4, FUNC, ...) FUNC
#define vec3_add2(x, y) _Generic((y), Vec3: vec3vec3_add, float: vec3float_add)(x, y)
#define vec3_add3(x, y, z) vec3_add2(vec3_add2(x, y), z)
#define vec3_add4(x, y, z, t) vec3_add2(vec3_add3(x, y, z), t)
#define vec3_add(...) _GET_MACRO(__VA_ARGS__, vec3_add4, vec3_add3, vec3_add2)(__VA_ARGS__)

#define vec3_sub(x, y) _Generic((y), Vec3: vec3vec3_sub, float: vec3float_sub)(x, y)
#define vec3_mul(x, y) _Generic((y), Vec3: vec3vec3_mul, float: vec3float_mul)(x, y)
#define vec3_div(x, y) _Generic((y), Vec3: vec3vec3_div, float: vec3float_div)(x, y)

struct Vec3 {
  float x;
  float y;
  float z;
};

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
Vec3 vec3_rand_between(float low, float high, PCG32State *rng);
Vec3 vec3_rand_unit_vector(PCG32State *rng);
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32State *rng);

struct Ray {
  Vec3 origin;
  Vec3 direction;
};

Vec3 ray_at(Ray ray, float t);

enum TextureType {
  SOLID_COLOR,
  CHECKER,
};

struct Texture {
  TextureType type;
  Vec3 color;
  float scale;
  Texture *even;
  Texture *odd;
};

Vec3 texture_value(Texture *texture, float u, float v, Vec3 p);

enum MaterialType {
  NORMAL,
  LAMBERTIAN,
  METAL,
  DIELECTRIC,
};

struct Material {
  MaterialType type;
  Texture albedo;
  float metal_fuzz;
  float eta;
};

struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material *material;
  float t;
  float u;
  float v;
  bool front_face;
};

struct Sphere {
  Vec3 center;
  float radius;
  Material *material;
};

struct World {
  Sphere *spheres;
  size_t n_spheres;
  Material *materials;
  size_t n_materials;
};

bool hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool hit_spheres(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);

struct PCG32State {
  uint64_t state;
  uint64_t inc;
};

void pcg32_srandom_r(PCG32State *rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_random_r(PCG32State *rng);
float pcg32_randomf_r(PCG32State *rng);

struct Camera {
  float aspect_ratio;
  int img_width;
  int img_height;
  int samples_per_pixel;
  int max_depth;
  float vfov;
  Vec3 look_from;
  Vec3 look_to;
  Vec3 vup;
  float dof_angle;
  float focal_length;
  Vec3 pixel00_loc;
  Vec3 pixel_delta_u;
  Vec3 pixel_delta_v;
  Vec3 u; // camera basis vectors
  Vec3 v;
  Vec3 w;
  Vec3 dof_disc_u;
  Vec3 dof_disc_v;
};

void camera_init(Camera *camera);
void camera_render(Camera *camera, World *world, uint8_t *buffer);
