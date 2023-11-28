#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>

typedef struct Material Material;

typedef struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material *material;
  float t;
  float u;
  float v;
  bool front_face;
} HitRecord;

typedef enum MaterialType {
  SURFACE_NORMAL,
  LAMBERTIAN,
  METAL,
  DIELECTRIC,
  DIFFUSE_LIGHT,
  ISOTROPIC,
} MaterialType;

struct Material {
  MaterialType tag;
  Texture *albedo;
  union {
    float fuzz; // for metal
    float eta;  // for dielectric
  };
};

bool Material_scatter(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, bool *skip_pdf, PCG32 *rng);
float Material_scatter_pdf(const Material *mat, Vec3 normal, Vec3 r_in, Vec3 r_out);
Vec3 Material_emit(const HitRecord *rec);

void SurfaceNormal_init(Material *self);
Material *SurfaceNormal_new();

void Lambertian_init(Material *self, Texture *albedo);
Material *Lambertian_new(Texture *albedo);

void Metal_init(Material *self, Texture *albedo, float fuzz);
Material *Metal_new(Texture *albedo, float fuzz);

void Dielectric_init(Material *self, float eta);
Material *Dielectric_new(float eta);

void DiffuseLight_init(Material *self, Texture *albedo);
Material *DiffuseLight_new(Texture *albedo);

void Isotropic_init(Material *self, Texture *albedo);
Material *Isotropic_new(Texture *albedo);

typedef struct ONB {
  Vec3 u;
  Vec3 v;
  Vec3 w;
} ONB;
Vec3 ONB_local(const ONB *self, Vec3 a);
void ONB_from_w(ONB *self, Vec3 w);

#endif // MATERIAL_H
