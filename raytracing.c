#include "raytracing.h"
#include <math.h>

Vec3 *vec3vec3_add_(Vec3 *self, Vec3 other) {
  self->x += other.x;
  self->y += other.y;
  self->z += other.z;
  return self;
}

Vec3 *vec3vec3_sub_(Vec3 *self, Vec3 other) {
  self->x -= other.x;
  self->y -= other.y;
  self->z -= other.z;
  return self;
}

Vec3 *vec3float_mul_(Vec3 *self, float other) {
  self->x *= other;
  self->y *= other;
  self->z *= other;
  return self;
}

Vec3 vec3_neg(Vec3 u) { return (Vec3){-u.x, -u.y, -u.z}; }
Vec3 vec3vec3_add(Vec3 u, Vec3 v) { return (Vec3){u.x + v.x, u.y + v.y, u.z + v.z}; }
Vec3 vec3vec3_sub(Vec3 u, Vec3 v) { return (Vec3){u.x - v.x, u.y - v.y, u.z - v.z}; }
Vec3 vec3vec3_mul(Vec3 u, Vec3 v) { return (Vec3){u.x * v.x, u.y * v.y, u.z * v.z}; }
Vec3 vec3float_add(Vec3 u, float v) { return (Vec3){u.x + v, u.y + v, u.z + v}; }
Vec3 vec3float_sub(Vec3 u, float v) { return (Vec3){u.x - v, u.y - v, u.z - v}; }
Vec3 vec3float_mul(Vec3 u, float v) { return (Vec3){u.x * v, u.y * v, u.z * v}; }

float vec3_length2(Vec3 u) { return vec3_dot(u, u); }
float vec3_length(Vec3 u) { return sqrtf(vec3_length2(u)); }
float vec3_dot(Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
Vec3 vec3_cross(Vec3 u, Vec3 v) { return (Vec3){u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y}; }
Vec3 vec3_unit(Vec3 u) { return vec3float_mul(u, 1.0f / vec3_length(u)); }
bool vec3_near_zero(Vec3 u) { return (fabsf(u.x) < 1e-8f) && (fabsf(u.y) < 1e-8f) && (fabsf(u.z) < 1e-8f); }

Vec3 ray_at(Ray ray, float t) { return vec3vec3_add(ray.origin, vec3float_mul(ray.direction, t)); }

bool hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  Vec3 oc = vec3_sub(ray->origin, sphere->center);
  float a = vec3_length2(ray->direction);
  float b = vec3_dot(oc, ray->direction);
  float c = vec3_length2(oc) - sphere->radius * sphere->radius;
  float disc = b * b - a * c;

  if (disc < 0)
    return false;

  float disc_sqrt = sqrtf(disc);
  float root = (-b - disc_sqrt) / a;
  if (root <= t_min || root >= t_max) {
    root = (-b + disc_sqrt) / a;
    if (root <= t_min || root >= t_max)
      return false;
  }

  hit_record->t = root;
  hit_record->p = ray_at(*ray, root);

  Vec3 outward_normal = vec3_mul(vec3_sub(hit_record->p, sphere->center), 1.0f / sphere->radius);
  hit_record->front_face = vec3_dot(ray->direction, outward_normal) < 0.0f;
  hit_record->normal = hit_record->front_face ? outward_normal : vec3_neg(outward_normal);
  hit_record->material = sphere->material;

  return true;
}

bool hit_spheres(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  bool hit_anything = false;

  for (int i = 0; i < world->n_spheres; i++)
    if (hit_sphere(world->spheres + i, ray, t_min, t_max, hit_record)) {
      t_max = hit_record->t;
      hit_anything = true;
    }

  return hit_anything;
}

// https://github.com/imneme/pcg-c-basic
void pcg32_srandom_r(PCG32State *rng, uint64_t initstate, uint64_t initseq) {
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg32_random_r(rng);
  rng->state += initstate;
  pcg32_random_r(rng);
}

uint32_t pcg32_random_r(PCG32State *rng) {
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + rng->inc;
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

float pcg32_randomf_r(PCG32State *rng) { return (pcg32_random_r(rng) >> 8) * (1.0f / 16777216.0f); }

Vec3 vec3_rand(PCG32State *rng) { return (Vec3){pcg32_randomf_r(rng), pcg32_randomf_r(rng), pcg32_randomf_r(rng)}; };
Vec3 vec3_rand_unit_vector(PCG32State *rng) {
  for (;;) {
    Vec3 v = vec3_rand(rng);
    float length2 = vec3_length2(v);
    if (length2 < 1)
      return vec3float_mul(v, 1.0f / sqrtf(length2));
  }
}
Vec3 vec3_rand_hemisphere(Vec3 normal, PCG32State *rng) {
  Vec3 v = vec3_rand_unit_vector(rng);
  return (vec3_dot(v, normal) > 0.0f) ? v : vec3_neg(v);
}

bool scatter_normal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = vec3float_mul(vec3float_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
  return false;
}

bool scatter_lambertian(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 new_direction = vec3vec3_add(hit_record->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(new_direction))
    new_direction = hit_record->normal;
  *scattered = new_direction;
  *color = hit_record->material->albedo;
  return true;
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3vec3_sub(incident, vec3float_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

bool scatter_metal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *scattered = vec3vec3_add(reflect(incident, hit_record->normal),
                            vec3float_mul(vec3_rand_unit_vector(rng), hit_record->material->metal_fuzz));
  *color = hit_record->material->albedo;
  return vec3_dot(*scattered, hit_record->normal) > 0; // check for degeneration
}

bool scatter_dielectric(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  float eta = hit_record->front_face ? 1.0f / hit_record->material->eta : hit_record->material->eta;
  incident = vec3_unit(incident);

  float cos_theta = fminf(-vec3_dot(incident, hit_record->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_randomf_r(rng)) {
    *scattered = reflect(incident, hit_record->normal);
  } else {
    Vec3 r_perp = vec3float_mul(vec3vec3_add(incident, vec3float_mul(hit_record->normal, cos_theta)), eta);
    Vec3 r_para = vec3float_mul(hit_record->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3vec3_add(r_perp, r_para);
  }
  *color = hit_record->material->albedo;
  return true;
}

bool scatter_default(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = (Vec3){0.0f, 0.0f, 0.0f};
  return false;
}

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  switch (hit_record->material->type) {
  case NORMAL:
    return scatter_normal(incident, hit_record, rng, scattered, color);
  case LAMBERTIAN:
    return scatter_lambertian(incident, hit_record, rng, scattered, color);
  case METAL:
    return scatter_metal(incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return scatter_dielectric(incident, hit_record, rng, scattered, color);
  default:
    return scatter_default(incident, hit_record, rng, scattered, color);
  }
}
