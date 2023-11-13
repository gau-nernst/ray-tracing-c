#include "material.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>

#define material_new(type, ...) define_struct_new(Material, type, __VA_ARGS__)

Material *SurfaceNormal_new() material_new(SURFACE_NORMAL);
Material *Lambertian_new(Texture albedo) material_new(LAMBERTIAN, albedo);
Material *Metal_new(Texture albedo, float fuzz) material_new(METAL, albedo, .fuzz = fuzz);
Material *Dielectric_new(Texture albedo, float eta) material_new(DIELECTRIC, albedo, .eta = eta);
Material *DiffuseLight_new(Texture albedo) material_new(DIFFUSE_LIGHT, albedo);

static bool lambertian_scatter(HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 new_direction = vec3_add(hit_record->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(new_direction))
    new_direction = hit_record->normal;
  *scattered = new_direction;
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

static Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

static bool metal_scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *scattered =
      vec3_add(reflect(incident, hit_record->normal), vec3_mul(vec3_rand_unit_vector(rng), hit_record->material->fuzz));
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return vec3_dot(*scattered, hit_record->normal) > 0; // check for degeneration
}

static bool dielectric_scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  float eta = hit_record->front_face ? 1.0f / hit_record->material->eta : hit_record->material->eta;
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
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  switch (hit_record->material->type) {
  case SURFACE_NORMAL:
    *color = vec3float_mul(vec3_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
    return false;
  case LAMBERTIAN:
    return lambertian_scatter(hit_record, rng, scattered, color);
  case METAL:
    return metal_scatter(incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return dielectric_scatter(incident, hit_record, rng, scattered, color);
  default:
    *color = (Vec3){0, 0, 0};
    return false;
  }
}

Vec3 emit(HitRecord *hit_record) {
  switch (hit_record->material->type) {
  case DIFFUSE_LIGHT:
    return *(Vec3 *)hit_record->material->albedo.ptr;
  default:
    return (Vec3){0, 0, 0};
  }
}
