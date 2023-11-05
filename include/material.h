#ifndef MATERIAL_H
#define MATERIAL_H

#include "vec3.h"

typedef struct Material Material;
typedef struct TextureChecker TextureChecker;
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

typedef struct Image {
  int width;
  int height;
  uint8_t *buffer;
} Image;

void image_load(Image *image, char *filename);

#define N_PERLIN 256

typedef struct PerlinNoise {
  float scale;
  Vec3 grad_field[N_PERLIN];
  int perm_x[N_PERLIN];
  int perm_y[N_PERLIN];
  int perm_z[N_PERLIN];
} PerlinNoise;

void perlin_noise_init(PerlinNoise *perlin, PCG32State *rng);

struct Texture {
  enum TextureType {
    SOLID,
    CHECKER,
    IMAGE,
    PERLIN,
  } type;
  union {
    Vec3 *color;
    TextureChecker *checker;
    Image *image;
    PerlinNoise *perlin;
  };
};

struct TextureChecker {
  float scale;
  Texture even;
  Texture odd;
};

Vec3 texture_value(Texture texture, float u, float v, Vec3 p);

struct Material {
  enum MaterialType {
    NORMAL,
    LAMBERTIAN,
    METAL,
    DIELECTRIC,
  } type;
  Texture albedo;
  float fuzz;
  float eta;
};

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);

#endif // MATERIAL_H
