#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "material.h"
#include "utils.h"
#include "vec3.h"
#include <stdbool.h>
#include <stdlib.h>

typedef struct Ray {
  Vec3 origin;
  Vec3 direction;
} Ray;

Vec3 ray_at(const Ray *ray, float t);

typedef enum HittableType {
  HITTABLE_LIST,
  SPHERE,
  QUAD,
} HittableType;

typedef struct Hittable {
  HittableType type;
  void *ptr;
} Hittable;

#define hittable(ptr)                                                                                                  \
  (Hittable) { _Generic((ptr), HittableList *: HITTABLE_LIST, Sphere *: SPHERE, Quad *: QUAD), ptr }

define_list_header(Hittable);

bool Hittable_hit(Hittable obj, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
bool HittableList_hit(const HittableList *list, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);

typedef struct Sphere {
  Vec3 center;
  float radius;
  Material material;
} Sphere;

Sphere *Sphere_new(Vec3 center, float radius, Material material);

typedef struct Quad {
  Vec3 Q;
  Vec3 u;
  Vec3 v;
  Vec3 normal;
  float D;
  Vec3 w;
  Material material;
} Quad;

void Quad_init(Quad *quad, Vec3 Q, Vec3 u, Vec3 v, Material mat);
Quad *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material mat);

HittableList *Box_new(Vec3 a, Vec3 b, Material mat);

// typedef struct AABB {
//   float x[3][2];
// } AABB;

typedef struct World {
  HittableList objects;
  MaterialList materials;
} World;

void World_init(World *world, size_t max_objects, size_t max_materials);

typedef struct Camera {
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
} Camera;

void Camera_init(Camera *camera);
void Camera_render(const Camera *camera, const World *world, uint8_t *buffer);

#endif // RAYTRACING_H
