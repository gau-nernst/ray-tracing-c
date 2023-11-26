#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "hittable.h"
#include "material.h"
#include "vec3.h"
#include <stddef.h>

typedef struct World {
  HittableList objects;
  Hittable *light;
} World;

void World_init(World *world, size_t max_objects);

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
