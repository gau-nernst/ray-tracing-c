#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "array.h"
#include "material.h"
#include "vec3.h"
#include <stdbool.h>
#include <stdlib.h>

typedef struct World World;
typedef struct Camera Camera;

typedef struct Ray {
  Vec3 origin;
  Vec3 direction;
} Ray;

Vec3 ray_at(const Ray *ray, float t);

typedef struct Sphere {
  Vec3 center;
  float radius;
  Material *material;
} Sphere;

Sphere sphere(Vec3 center, float radius, Material *material);

typedef struct Quad {
  Vec3 Q;
  Vec3 u;
  Vec3 v;
  Vec3 normal;
  float D;
  Vec3 w;
  Material *material;
} Quad;

void quad_init(Quad *quad, Vec3 Q, Vec3 u, Vec3 v, Material *material);

typedef struct AABB {
  float x[3][2];
} AABB;

define_array_header(Vec3);
define_array_header(Sphere);
define_array_header(Quad);
define_array_header(Material);
define_array_header(Checker);
define_array_header(Image);
define_array_header(Perlin);

#define array_append(array, obj)                                                                                       \
  _Generic((array),                                                                                                    \
      Vec3Array *: Vec3Array_append,                                                                                   \
      SphereArray *: SphereArray_append,                                                                               \
      MaterialArray *: MaterialArray_append,                                                                           \
      CheckerArray *: CheckerArray_append,                                                                             \
      ImageArray *: ImageArray_append,                                                                                 \
      PerlinArray *: PerlinArray_append)(array, obj)

#define array_next(array)                                                                                              \
  _Generic((array),                                                                                                    \
      Vec3Array *: Vec3Array_next,                                                                                     \
      SphereArray *: SphereArray_next,                                                                                 \
      MaterialArray *: MaterialArray_next,                                                                             \
      CheckerArray *: CheckerArray_next,                                                                               \
      ImageArray *: ImageArray_next,                                                                                   \
      PerlinArray *: PerlinArray_next)(array)

struct World {
  SphereArray spheres;
  QuadArray quads;
  Vec3Array colors;
  MaterialArray materials;
  CheckerArray checkers;
  ImageArray images;
  PerlinArray perlins;
};

bool hit_objects(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);

struct Camera {
  float aspect_ratio;
  int img_width;
  int img_height;
  int samples_per_pixel;
  int max_depth;
  Vec3 background;
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
void camera_render(const Camera *camera, const World *world, uint8_t *buffer);

#endif // RAYTRACING_H
