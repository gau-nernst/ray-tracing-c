#include "vec3.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#define clamp(x, lo, hi) min(max(x, lo), hi)

typedef struct Ray Ray;
typedef struct HitRecord HitRecord;
typedef struct Sphere Sphere;
typedef enum MaterialType MaterialType;
typedef struct Material Material;
typedef struct World World;
typedef struct Camera Camera;
typedef enum TextureType TextureType;
typedef struct Texture Texture;

struct Ray {
  Vec3 origin;
  Vec3 direction;
};

Vec3 ray_at(Ray ray, float t);

enum TextureType {
  SOLID_COLOR,
  CHECKER,
};

struct Texture {
  TextureType type;
  Vec3 color;
  float scale;
  Texture *even;
  Texture *odd;
};

Vec3 texture_value(Texture *texture, float u, float v, Vec3 p);

enum MaterialType {
  NORMAL,
  LAMBERTIAN,
  METAL,
  DIELECTRIC,
};

struct Material {
  MaterialType type;
  Texture *albedo;
  float metal_fuzz;
  float eta;
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

struct Sphere {
  Vec3 center;
  float radius;
  Material *material;
};

struct World {
  Sphere *spheres;
  size_t n_spheres;
  Material *materials;
  size_t n_materials;
  Texture *textures;
  size_t n_textures;
};

bool hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool hit_spheres(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color);

struct Camera {
  float aspect_ratio;
  int img_width;
  int img_height;
  int samples_per_pixel;
  int max_depth;
  float vfov;
  Vec3 look_from;
  Vec3 look_to;
  Vec3 vup;
  float dof_angle;
  float focal_length;
  Vec3 pixel00_loc;
  Vec3 pixel_delta_u;
  Vec3 pixel_delta_v;
  Vec3 u; // camera basis vectors
  Vec3 v;
  Vec3 w;
  Vec3 dof_disc_u;
  Vec3 dof_disc_v;
};

void camera_init(Camera *camera);
void camera_render(Camera *camera, World *world, uint8_t *buffer);
