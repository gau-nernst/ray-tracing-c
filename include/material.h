#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>

typedef struct Material Material;
typedef struct HitRecord HitRecord;

struct Material {
  bool (*scatter)(HitRecord *rec, Vec3 incident, PCG32 *rng, Vec3 *scattered, Vec3 *color);
  Vec3 (*emit)(HitRecord *rec);
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

typedef struct SurfaceNormal {
  Material mat;
} SurfaceNormal;
void SurfaceNormal_init(SurfaceNormal *self);
Material *SurfaceNormal_new();

typedef struct Lambertian {
  Material mat;
  Texture *albedo;
} Lambertian;
void Lambertian_init(Lambertian *self, Texture *albedo);
Material *Lambertian_new(Texture *albedo);

typedef struct Metal {
  Material mat;
  Texture *albedo;
  float fuzz;
} Metal;
void Metal_init(Metal *self, Texture *albedo, float fuzz);
Material *Metal_new(Texture *albedo, float fuzz);

typedef struct Dielectric {
  Material mat;
  float eta;
} Dielectric;
void Dielectric_init(Dielectric *self, float eta);
Material *Dielectric_new(float eta);

typedef struct DiffuseLight {
  Material mat;
  Texture *albedo;
} DiffuseLight;
void DiffuseLight_init(DiffuseLight *self, Texture *albedo);
Material *DiffuseLight_new(Texture *albedo);

typedef struct Isotropic {
  Material mat;
  Texture *albedo;
} Isotropic;
void Isotropic_init(Isotropic *self, Texture *albedo);
Material *Isotropic_new(Texture *albedo);

#endif // MATERIAL_H
