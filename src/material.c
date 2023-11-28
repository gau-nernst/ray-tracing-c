#include "material.h"
#include "utils.h"
#include "vec3.h"
#include <math.h>
#include <stdlib.h>

static Vec3 _Texture_value(const HitRecord *rec) {
  const Texture *texture = rec->material->albedo;
  return texture->value(texture, rec->u, rec->v, rec->p);
}

#define define_material_new(Type, ...)                                                                                 \
  {                                                                                                                    \
    Material *obj = my_malloc(sizeof(Material));                                                                       \
    Type##_init(obj, ##__VA_ARGS__);                                                                                   \
    return obj;                                                                                                        \
  }

void SurfaceNormal_init(Material *self) { self->tag = SURFACE_NORMAL; }
Material *SurfaceNormal_new() define_material_new(SurfaceNormal);

// sample from p(theta) = cos(theta) / pi, 0 <= theta <= pi/2
static Vec3 rand_cosine_theta(PCG32 *rng) {
  float r1 = pcg32_f32(rng);
  float r2 = pcg32_f32(rng);
  float phi = 2.0f * (float)M_PI * r1;
  return vec3(cosf(phi) * sqrtf(r2), sinf(phi) * sqrtf(r2), sqrtf(1.0f - r2));
}
static bool Lambertian_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng) {
  ONB onb;
  ONB_from_w(&onb, rec->normal);

  *r_out = ONB_local(&onb, rand_cosine_theta(rng));
  *color = _Texture_value(rec);
  *skip_pdf = false;
  return true;
}
static float Lambertian_scatter_pdf(Vec3 normal, Vec3 r_in, Vec3 r_out) {
  float cos_theta = vec3_dot(normal, vec3_normalize(r_out));
  return cos_theta < 0.0f ? 0.0f : cos_theta / (float)M_PI;
}
void Lambertian_init(Material *self, Texture *albedo) { *self = (Material){LAMBERTIAN, albedo}; }
Material *Lambertian_new(Texture *albedo) define_material_new(Lambertian, albedo);

static Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}
static bool Metal_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng) {
  Vec3 reflected = reflect(vec3_normalize(r_in), rec->normal);
  *r_out = vec3_add(reflected, vec3_mul(vec3_rand_unit_vector(rng), rec->material->fuzz));
  *color = _Texture_value(rec);
  *skip_pdf = true;
  return vec3_dot(*r_out, rec->normal) > 0.0f; // check for degeneration
}
void Metal_init(Material *self, Texture *albedo, float fuzz) { *self = (Material){METAL, albedo, .fuzz = fuzz}; }
Material *Metal_new(Texture *albedo, float fuzz) define_material_new(Metal, albedo, fuzz);

static bool Dielectric_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng) {
  float eta = rec->material->eta;
  if (rec->front_face)
    eta = 1.0f / eta;
  r_in = vec3_normalize(r_in);

  float cos_theta = fminf(-vec3_dot(r_in, rec->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_f32(rng)) {
    *r_out = reflect(r_in, rec->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(r_in, vec3_mul(rec->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(rec->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *r_out = vec3_add(r_perp, r_para);
  }
  *color = vec3(1, 1, 1);
  *skip_pdf = true;
  return true;
}
void Dielectric_init(Material *self, float eta) { *self = (Material){DIELECTRIC, NULL, .eta = eta}; }
Material *Dielectric_new(float eta) define_material_new(Dielectric, eta);

void DiffuseLight_init(Material *self, Texture *albedo) { *self = (Material){DIFFUSE_LIGHT, albedo}; }
Material *DiffuseLight_new(Texture *albedo) define_material_new(DiffuseLight, albedo);

static bool Isotropic_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng) {
  *r_out = vec3_rand_unit_vector(rng);
  *color = _Texture_value(rec);
  *skip_pdf = false;
  return true;
}
static float Isotropic_scatter_pdf(Vec3 normal, Vec3 r_in, Vec3 r_out) { return 1.0f / (4.0f * (float)M_PI); }
void Isotropic_init(Material *self, Texture *albedo) { *self = (Material){ISOTROPIC, albedo}; }
Material *Isotropic_new(Texture *albedo) define_material_new(Isotropic, albedo);

bool Material_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng) {
  switch (rec->material->tag) {
  case SURFACE_NORMAL:
    *color = vec3float_mul(vec3_add(rec->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
    *skip_pdf = true;
    return false;
  case LAMBERTIAN:
    return Lambertian_scatter(rec, r_in, r_out, color, skip_pdf, rng);
  case METAL:
    return Metal_scatter(rec, r_in, r_out, color, skip_pdf, rng);
  case DIELECTRIC:
    return Dielectric_scatter(rec, r_in, r_out, color, skip_pdf, rng);
  case DIFFUSE_LIGHT:
    *skip_pdf = true;
    return false;
  case ISOTROPIC:
    return Isotropic_scatter(rec, r_in, r_out, color, skip_pdf, rng);
  }
}

float Material_scatter_pdf(const Material *mat, Vec3 normal, Vec3 r_in, Vec3 r_out) {
  switch (mat->tag) {
  case LAMBERTIAN:
    return Lambertian_scatter_pdf(normal, r_in, r_out);
  case ISOTROPIC:
    return Isotropic_scatter_pdf(normal, r_in, r_out);
  default:
    return 0.0f;
  }
}

Vec3 Material_emit(const HitRecord *rec) {
  if (rec->material->tag == DIFFUSE_LIGHT)
    return rec->front_face ? _Texture_value(rec) : VEC3_ZERO;
  return VEC3_ZERO;
}

Vec3 ONB_local(const ONB *self, Vec3 a) {
  return vec3_add(vec3_mul(self->u, a.x), vec3_mul(self->v, a.y), vec3_mul(self->w, a.z));
}
void ONB_from_w(ONB *self, Vec3 w) {
  self->w = vec3_normalize(w);
  Vec3 a = fabsf(self->w.x) > 0.9f ? vec3(0, 1, 0) : vec3(1, 0, 0);
  self->v = vec3_normalize(vec3_cross(self->w, a));
  self->u = vec3_cross(self->w, self->v);
}
