#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "hittable.h"
#include "material.h"
#include "vec3.h"
#include <stddef.h>

typedef struct World {
  HittableList objects;
  HittableList lights;
} World;

void World_init(World *world, size_t max_objects);

typedef struct Camera {
  Vec3 background;
  float vfov; // may remove?
  Vec3 look_from;
  Vec3 look_to;
  Vec3 vup; // may remove?
  float dof_angle;
  float focal_length;
  Vec3 u; // camera basis vectors
  Vec3 v;
  Vec3 w;
  Vec3 dof_disc_u;
  Vec3 dof_disc_v;
} Camera;

void Camera_init(Camera *self, Vec3 background, float vfov, float focal_length, Vec3 look_from, Vec3 look_to,
                 float dof_angle);

typedef struct Renderer {
  int img_width;
  int img_height;
  int samples_per_pixel;
  int max_depth;
  float lights_sampling_prob;
} Renderer;

void Renderer_init(Renderer *self, int img_width, float aspect_ratio, int samples_per_pixel, int max_depth,
                   float lights_sampling_prob);

void render(const Renderer *renderer, const Camera *camera, const World *world, uint8_t *buffer);

#endif // RAYTRACING_H
