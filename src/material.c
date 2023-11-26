#include "material.h"
#include "utils.h"
#include "vec3.h"
#include <math.h>
#include <stdlib.h>

static Vec3 _Texture_value(Texture *texture, HitRecord *rec) { return texture->value(texture, rec->u, rec->v, rec->p); }
static Vec3 Material_emit(HitRecord *rec) { return VEC3_ZERO; }

static bool SurfaceNormal_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  *color = vec3float_mul(vec3_add(rec->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
  return false;
}
void SurfaceNormal_init(SurfaceNormal *self) { *self = (SurfaceNormal){SurfaceNormal_scatter, Material_emit}; }
Material *SurfaceNormal_new() define_init_new(SurfaceNormal);

static bool Lambertian_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  Lambertian *mat = (Lambertian *)rec->material;
  *scattered = vec3_add(rec->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(*scattered)) // remove degenerate rays
    *scattered = rec->normal;
  *color = _Texture_value(mat->albedo, rec);
  return true;
}
void Lambertian_init(Lambertian *self, Texture *albedo) {
  *self = (Lambertian){{Lambertian_scatter, Material_emit}, albedo};
}
Material *Lambertian_new(Texture *albedo) define_init_new(Lambertian, albedo);

static Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}
static bool Metal_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  Metal *mat = (Metal *)rec->material;
  Vec3 reflected = reflect(vec3_normalize(incident), rec->normal);
  *scattered = vec3_add(reflected, vec3_mul(vec3_rand_unit_vector(rng), mat->fuzz));
  *color = _Texture_value(mat->albedo, rec);
  return vec3_dot(*scattered, rec->normal) > 0.0f; // check for degeneration
}
void Metal_init(Metal *self, Texture *albedo, float fuzz) {
  *self = (Metal){{Metal_scatter, Material_emit}, albedo, fuzz};
}
Material *Metal_new(Texture *albedo, float fuzz) define_init_new(Metal, albedo, fuzz);

static bool Dielectric_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  Dielectric *mat = (Dielectric *)rec->material;
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
  *color = vec3(1, 1, 1);
  return true;
}
void Dielectric_init(Dielectric *self, float eta) { *self = (Dielectric){{Dielectric_scatter, Material_emit}, eta}; }
Material *Dielectric_new(float eta) define_init_new(Dielectric, eta);

static bool DiffuseLight_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  return false;
}
static Vec3 DiffuseLight_emit(HitRecord *rec) {
  DiffuseLight *mat = (DiffuseLight *)rec->material;
  return _Texture_value(mat->albedo, rec);
}
void DiffuseLight_init(DiffuseLight *self, Texture *albedo) {
  *self = (DiffuseLight){{DiffuseLight_scatter, DiffuseLight_emit}, albedo};
}
Material *DiffuseLight_new(Texture *albedo) define_init_new(DiffuseLight, albedo);

static bool Isotropic_scatter(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color) {
  Isotropic *mat = (Isotropic *)rec->material;
  *scattered = vec3_rand_unit_vector(rng);
  *color = _Texture_value(mat->albedo, rec);
  return true;
}
void Isotropic_init(Isotropic *self, Texture *albedo) {
  *self = (Isotropic){{Isotropic_scatter, Material_emit}, albedo};
}
Material *Isotropic_new(Texture *albedo) define_init_new(Isotropic, albedo);
