#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "material.h"
#include "vec3.h"
#include <stdbool.h>
#include <stdlib.h>

typedef struct Ray Ray;
typedef struct Sphere Sphere;
typedef struct World World;
typedef struct Camera Camera;

struct Ray {
  Vec3 origin;
  Vec3 direction;
};

Vec3 ray_at(const Ray *ray, float t);

struct Sphere {
  Vec3 center;
  float radius;
  Material *material;
};

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

struct World {
  size_t n_spheres;
  size_t n_quads;
  size_t n_materials;
  size_t n_colors;
  size_t n_checkers;
  size_t n_images;
  size_t n_perlins;

  Sphere *spheres;
  Quad *quads;
  Material *materials;
  Vec3 *colors;
  Checker *checkers;
  Image *images;
  Perlin *perlins;
};

void world_malloc(World *world);
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
