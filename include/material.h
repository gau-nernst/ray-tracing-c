#ifndef MATERIAL_H
#define MATERIAL_H

#include "vec3.h"

typedef struct Material Material;
typedef struct Checker Checker;
typedef struct Texture Texture;

typedef struct HitRecord {
  Vec3 p;
  Vec3 normal;
  Material *material;
  float t;
  float u;
  float v;
  bool front_face;
} HitRecord;

typedef struct Image {
  int width;
  int height;
  uint8_t *buffer;
} Image;

void image_load(Image *image, char *filename);

#define N_PERLIN 256

typedef struct Perlin {
  float scale;
  int depth;
  Vec3 grad_field[N_PERLIN];
  int perm_x[N_PERLIN];
  int perm_y[N_PERLIN];
  int perm_z[N_PERLIN];
} Perlin;

void perlin_init(Perlin *perlin, PCG32State *rng);

struct Texture {
  enum TextureType {
    SOLID,
    CHECKER,
    IMAGE,
    PERLIN,
  } type;
  void *ptr;
};

Texture texture(enum TextureType type, void *ptr);

struct Checker {
  float scale;
  Texture even;
  Texture odd;
};

Vec3 texture_value(Texture texture, float u, float v, Vec3 p);

struct Material {
  enum MaterialType {
    SURFACE_NORMAL,
    LAMBERTIAN,
    METAL,
    DIELECTRIC,
    DIFFUSE_LIGHT,
  } type;
  Texture albedo;
  union {
    float fuzz; // used by METAL
    float eta;  // used by DIELECTRIC
  };
};

Material surface_normal();

#define lambertian(albedo) _Generic((albedo), Vec3 *: lambertian_color, Texture: lambertian_texture)(albedo)
Material lambertian_color(Vec3 *color);
Material lambertian_texture(Texture albedo);

#define metal(albedo, fuzz) _Generic((albedo), Vec3 *: metal_color, Texture: metal_texture)(albedo, fuzz)
Material metal_color(Vec3 *color, float fuzz);
Material metal_texture(Texture albedo, float fuzz);

#define dielectric(albedo, eta) _Generic((albedo), Vec3 *: dielectric_color, Texture: dielectric_texture)(albedo, eta)
Material dielectric_color(Vec3 *color, float eta);
Material dielectric_texture(Texture albedo, float eta);

#define diffuse_light(albedo) _Generic((albedo), Vec3 *: diffuse_light_color, Texture: diffuse_light_texture)(albedo)
Material diffuse_light_color(Vec3 *color);
Material diffuse_light_texture(Texture albedo);

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);
Vec3 emit(HitRecord *hit_record);

#endif // MATERIAL_H
