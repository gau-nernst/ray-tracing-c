#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>

typedef enum MaterialType {
  SURFACE_NORMAL,
  LAMBERTIAN,
  METAL,
  DIELECTRIC,
  DIFFUSE_LIGHT,
  ISOTROPIC,
} MaterialType;

typedef struct Material {
  MaterialType type;
  void *ptr;
} Material;

define_list_header(Material);

#define material(ptr)                                                                                                  \
  (Material) {                                                                                                         \
    _Generic((ptr),                                                                                                    \
        SurfaceNormal *: SURFACE_NORMAL,                                                                               \
        Lambertian *: LAMBERTIAN,                                                                                      \
        Metal *: METAL,                                                                                                \
        Dielectric *: DIELECTRIC,                                                                                      \
        DiffuseLight *: DIFFUSE_LIGHT,                                                                                 \
        Isotropic *: ISOTROPIC),                                                                                       \
        ptr                                                                                                            \
  }

typedef struct SurfaceNormal SurfaceNormal;

typedef struct Lambertian {
  Texture *albedo;
} Lambertian;

typedef struct Metal {
  Texture *albedo;
  float fuzz;
} Metal;

typedef struct Dielectric {
  Texture *albedo;
  float eta;
} Dielectric;

typedef struct DiffuseLight {
  Texture *albedo;
} DiffuseLight;

typedef struct Isotropic {
  Texture *albedo;
} Isotropic;

SurfaceNormal *SurfaceNormal_new();
Lambertian *Lambertian_new(Texture *albedo);
Metal *Metal_new(Texture *albedo, float fuzz);
Dielectric *Dielectric_new(Texture *albedo, float eta);
DiffuseLight *DiffuseLight_new(Texture *albedo);
Isotropic *Isotropic_new(Texture *albedo);

typedef struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material material;
  float t;
  float u;
  float v;
  bool front_face;
} HitRecord;

bool scatter(Vec3 incident, HitRecord *rec, PCG32State *rng, Vec3 *scattered, Vec3 *color);
Vec3 emit(HitRecord *rec);

#endif // MATERIAL_H
