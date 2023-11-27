#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>

typedef struct Material Material;
typedef struct HitRecord HitRecord;
typedef struct Ray Ray;

typedef struct MaterialVTable {
  bool (*scatter)(const HitRecord *rec, Vec3 r_in, Vec3 *r_out, Vec3 *color, float *pdf, PCG32 *rng);
  float (*scatter_pdf)(const HitRecord *rec, Vec3 r_in, Vec3 r_out);
  Vec3 (*emit)(const HitRecord *rec);
} MaterialVTable;

struct Material {
  MaterialVTable *vtable;
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
Vec3 ONB_local(const ONB *self, Vec3 a);
void ONB_from_w(ONB *self, Vec3 w);

#endif // MATERIAL_H
