#include "material.h"
#include <math.h>
#include <stdlib.h>

define_list_source(Material);

#define material_new(type, ...) define_struct_new(Material, type, __VA_ARGS__)

SurfaceNormal *SurfaceNormal_new() { return NULL; };
Lambertian *Lambertian_new(Texture albedo) define_struct_new(Lambertian, albedo);
Metal *Metal_new(Texture albedo, float fuzz) define_struct_new(Metal, albedo, fuzz);
Dielectric *Dielectric_new(Texture albedo, float eta) define_struct_new(Dielectric, albedo, eta);
DiffuseLight *DiffuseLight_new(Texture albedo) define_struct_new(DiffuseLight, albedo);
Isotropic *Isotropic_new(Texture albedo) define_struct_new(Isotropic, albedo);

static bool lambertian_scatter(Lambertian *mat, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);
static bool metal_scatter(Metal *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                          Vec3 *color);
static bool dielectric_scatter(Dielectric *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                               Vec3 *color);
static bool isotropic_scatter(Isotropic *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                              Vec3 *color);

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Material mat = hit_record->material;
  switch (mat.type) {
  case SURFACE_NORMAL:
    *color = vec3float_mul(vec3_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
    return false;
  case LAMBERTIAN:
    return lambertian_scatter(mat.ptr, hit_record, rng, scattered, color);
  case METAL:
    return metal_scatter(mat.ptr, incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return dielectric_scatter(mat.ptr, incident, hit_record, rng, scattered, color);
  case ISOTROPIC:
    return isotropic_scatter(mat.ptr, incident, hit_record, rng, scattered, color);
  default:
    *color = (Vec3){0, 0, 0};
    return false;
  }
}

static bool lambertian_scatter(Lambertian *mat, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 new_direction = vec3_add(hit_record->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(new_direction))
    new_direction = hit_record->normal;
  *scattered = new_direction;
  *color = texture_value(mat->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

static Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

static bool metal_scatter(Metal *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                          Vec3 *color) {
  *scattered = vec3_add(reflect(incident, hit_record->normal), vec3_mul(vec3_rand_unit_vector(rng), mat->fuzz));
  *color = texture_value(mat->albedo, hit_record->u, hit_record->v, hit_record->p);
  return vec3_dot(*scattered, hit_record->normal) > 0; // check for degeneration
}

static bool dielectric_scatter(Dielectric *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                               Vec3 *color) {
  float eta = hit_record->front_face ? 1.0f / mat->eta : mat->eta;
  incident = vec3_unit(incident);

  float cos_theta = fminf(-vec3_dot(incident, hit_record->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_f32(rng)) {
    *scattered = reflect(incident, hit_record->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(incident, vec3_mul(hit_record->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(hit_record->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3_add(r_perp, r_para);
  }
  *color = texture_value(mat->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

static bool isotropic_scatter(Isotropic *mat, Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered,
                              Vec3 *color) {
  *scattered = vec3_rand_unit_vector(rng);
  *color = texture_value(mat->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

static Vec3 DiffuseLight_emit(DiffuseLight *mat, HitRecord *hit_record) {
  return texture_value(mat->albedo, hit_record->u, hit_record->v, hit_record->p);
}

Vec3 emit(HitRecord *hit_record) {
  Material mat = hit_record->material;
  switch (mat.type) {
  case DIFFUSE_LIGHT:
    return DiffuseLight_emit(mat.ptr, hit_record);
  default:
    return (Vec3){0, 0, 0};
  }
}
