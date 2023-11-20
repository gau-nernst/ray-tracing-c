#include "material.h"
#include <math.h>
#include <stdlib.h>

define_list_source(Material);

#define material_new(type, ...) define_struct_new(Material, type, __VA_ARGS__)

SurfaceNormal *SurfaceNormal_new() { return NULL; };
Lambertian *Lambertian_new(Texture *albedo) define_struct_new(Lambertian, albedo);
Metal *Metal_new(Texture *albedo, float fuzz) define_struct_new(Metal, albedo, fuzz);
Dielectric *Dielectric_new(Texture *albedo, float eta) define_struct_new(Dielectric, albedo, eta);
DiffuseLight *DiffuseLight_new(Texture *albedo) define_struct_new(DiffuseLight, albedo);
Isotropic *Isotropic_new(Texture *albedo) define_struct_new(Isotropic, albedo);

static bool Lambertian_scatter(Lambertian *mat, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color);
static bool Metal_scatter(Metal *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color);
static bool Dielectric_scatter(Dielectric *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered,
                               Vec3 *color);
static bool Isotropic_scatter(Isotropic *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered,
                              Vec3 *color);

bool scatter(Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Material mat = rec->material;
  switch (mat.type) {
  case SURFACE_NORMAL:
    *color = vec3float_mul(vec3_add(rec->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
    return false;
  case LAMBERTIAN:
    return Lambertian_scatter(mat.ptr, rec, rng, scattered, color);
  case METAL:
    return Metal_scatter(mat.ptr, incident, rec, rng, scattered, color);
  case DIELECTRIC:
    return Dielectric_scatter(mat.ptr, incident, rec, rng, scattered, color);
  case ISOTROPIC:
    return Isotropic_scatter(mat.ptr, incident, rec, rng, scattered, color);
  default:
    *color = VEC3_ZERO;
    return false;
  }
}

static Vec3 _Texture_value(Texture *texture, HitRecord *rec) { return texture->value(texture, rec->u, rec->v, rec->p); }

static bool Lambertian_scatter(Lambertian *mat, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *scattered = vec3_add(rec->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(*scattered)) // remove degenerate rays
    *scattered = rec->normal;
  *color = _Texture_value(mat->albedo, rec);
  return true;
}

static Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

static bool Metal_scatter(Metal *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 reflected = reflect(vec3_normalize(incident), rec->normal);
  *scattered = vec3_add(reflected, vec3_mul(vec3_rand_unit_vector(rng), mat->fuzz));
  *color = _Texture_value(mat->albedo, rec);
  return vec3_dot(*scattered, rec->normal) > 0.0f; // check for degeneration
}

static bool Dielectric_scatter(Dielectric *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered,
                               Vec3 *color) {
  float eta = rec->front_face ? 1.0f / mat->eta : mat->eta;
  incident = vec3_normalize(incident);

  float cos_theta = fminf(-vec3_dot(incident, rec->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_f32(rng)) {
    *scattered = reflect(incident, rec->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(incident, vec3_mul(rec->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(rec->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3_add(r_perp, r_para);
  }
  *color = _Texture_value(mat->albedo, rec);
  return true;
}

static bool Isotropic_scatter(Isotropic *mat, Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered,
                              Vec3 *color) {
  *scattered = vec3_rand_unit_vector(rng);
  *color = _Texture_value(mat->albedo, rec);
  return true;
}

static Vec3 DiffuseLight_emit(DiffuseLight *mat, HitRecord *rec) { return _Texture_value(mat->albedo, rec); }

Vec3 emit(HitRecord *rec) {
  Material mat = rec->material;
  switch (mat.type) {
  case DIFFUSE_LIGHT:
    return DiffuseLight_emit(mat.ptr, rec);
  default:
    return VEC3_ZERO;
  }
}
