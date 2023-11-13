#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "vec3.h"

typedef struct Material {
  enum MaterialType {
    SURFACE_NORMAL,
    LAMBERTIAN,
    METAL,
    DIELECTRIC,
    DIFFUSE_LIGHT,
  } type;
  Texture albedo;
  union {
    float fuzz; // for metal
    float eta;  // for dielectric
  };
} Material;

Material *SurfaceNormal_new();
Material *Lambertian_new(Texture albedo);
Material *Metal_new(Texture albedo, float fuzz);
Material *Dielectric_new(Texture albedo, float eta);
Material *DiffuseLight_new(Texture albedo);

typedef struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material *material;
  float t;
  float u;
  float v;
  bool front_face;
} HitRecord;

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);
Vec3 emit(HitRecord *hit_record);

#endif // MATERIAL_H
