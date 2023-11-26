#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>

typedef struct Material Material;
typedef struct HitRecord HitRecord;
typedef struct Ray Ray;

struct Material {
  bool (*scatter)(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, PCG32 *rng);
  Vec3 (*emit)(const HitRecord *rec);
  float (*scattering_pdf)(const HitRecord *rec, Vec3 r_in, Vec3 r_out);
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

typedef struct ONB {
  Vec3 u;
  Vec3 v;
  Vec3 w;
} ONB;
Vec3 ONB_local(const ONB *self, float u, float v, float w);
void ONB_from_w(ONB *self, Vec3 w);

#endif // MATERIAL_H
