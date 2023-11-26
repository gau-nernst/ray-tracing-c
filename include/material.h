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
  Texture *albedo;
  union {
    float fuzz; // for metal
    float eta;  // for dielectric
  };
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

#endif // MATERIAL_H
