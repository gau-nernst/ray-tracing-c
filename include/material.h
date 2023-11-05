#ifndef MATERIAL_H
#define MATERIAL_H

#include "vec3.h"

typedef struct Material Material;
typedef struct Texture Texture;
typedef struct HitRecord HitRecord;

struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material *material;
  float t;
  float u;
  float v;
  bool front_face;
};

struct Texture {
  enum TextureType {
    SOLID,
    CHECKER,
    IMAGE,
    PERLIN,
  } type;
  union {
    Vec3 color; // SOLID
    struct {    // CHECKER
      float checker_scale;
      Texture *even;
      Texture *odd;
    };
    struct { // IMAGE
      int width;
      int height;
      uint8_t *image;
    };
    struct { // PERLIN
      float perlin_scale;
      float *values;
      int *perm_x;
      int *perm_y;
      int *perm_z;
    };
  };
};

void texture_image_load(Texture *texture, char *filename);
void texture_perlin_init(Texture *texture, PCG32State *rng);
Vec3 texture_value(Texture *texture, float u, float v, Vec3 p);

struct Material {
  enum MaterialType {
    NORMAL,
    LAMBERTIAN,
    METAL,
    DIELECTRIC,
  } type;
  Texture *albedo;
  float fuzz;
  float eta;
};

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);

#endif // MATERIAL_H
