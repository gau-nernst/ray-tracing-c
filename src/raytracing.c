#include "raytracing.h"
#define _USE_MATH_DEFINES // for MSVC
#include <math.h>
#include <stdio.h>

Vec3 ray_at(Ray ray, float t) { return vec3vec3_add(ray.origin, vec3float_mul(ray.direction, t)); }

bool sphere_hit(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  Vec3 oc = vec3_sub(ray->origin, sphere->center);
  float a = vec3_length2(ray->direction);
  float b = vec3_dot(oc, ray->direction);
  float c = vec3_length2(oc) - sphere->radius * sphere->radius;
  float disc = b * b - a * c;

  if (disc < 0)
    return false;

  float disc_sqrt = sqrtf(disc);
  float root = (-b - disc_sqrt) / a;
  if (root <= t_min || root >= t_max) {
    root = (-b + disc_sqrt) / a;
    if (root <= t_min || root >= t_max)
      return false;
  }

  hit_record->t = root;
  hit_record->p = ray_at(*ray, root);

  Vec3 outward_normal = vec3_div(vec3_sub(hit_record->p, sphere->center), sphere->radius);
  hit_record->front_face = vec3_dot(ray->direction, outward_normal) < 0.0f;
  hit_record->normal = hit_record->front_face ? outward_normal : vec3_neg(outward_normal);
  hit_record->u = (atan2f(-outward_normal.x[2], outward_normal.x[0]) + M_PI) * M_1_PI * 0.5;
  hit_record->v = acosf(-outward_normal.x[1]) * M_1_PI;
  hit_record->material = sphere->material;

  return true;
}

void quad_init(Quad *quad) {
  quad->normal = vec3_unit(vec3_cross(quad->u, quad->v));
  quad->D = vec3_dot(quad->normal, quad->Q);
}

bool quad_hit(const Quad *quad, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) { return false; }

bool aabb_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max) {
  for (int a = 0; a < 3; a++) {
    float invD = 1.0f / ray->direction.x[a];
    float t0 = (min(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;
    float t1 = (max(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;

    if (invD < 0) {
      float tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    if (t0 > t_min)
      t_min = t0;
    if (t1 < t_max)
      t_max = t1;
    if (t_max < t_min)
      return false;
  }
  return true;
}

AABB aabb_pad(const AABB *aabb) {
  float delta = 1e-4f;
  AABB padded;
  for (int a = 0; a < 3; a++) {
    if (aabb->x[a][1] - aabb->x[a][0] < delta) {
      padded.x[a][0] = aabb->x[a][0] - delta;
      padded.x[a][1] = aabb->x[a][1] + delta;
    } else {
      padded.x[a][0] = aabb->x[a][0];
      padded.x[a][1] = aabb->x[a][1];
    }
  }
  return padded;
}

void world_init(World *world) {
  try_malloc(world->spheres, sizeof(Sphere) * world->n_spheres);
  try_malloc(world->quads, sizeof(Quad) * world->n_quads);
  try_malloc(world->materials, sizeof(Material) * world->n_materials);
  try_malloc(world->colors, sizeof(Vec3) * world->n_colors);
  try_malloc(world->checkers, sizeof(Checker) * world->n_checkers);
  try_malloc(world->images, sizeof(Image) * world->n_images);
  try_malloc(world->perlins, sizeof(Perlin) * world->n_perlins);
}

bool hit_objects(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  bool hit_anything = false;

  for (int i = 0; i < world->n_spheres; i++)
    if (sphere_hit(world->spheres + i, ray, t_min, t_max, hit_record)) {
      t_max = hit_record->t;
      hit_anything = true;
    }

  for (int i = 0; i < world->n_quads; i++)
    if (quad_hit(world->quads + i, ray, t_min, t_max, hit_record)) {
      t_max = hit_record->t;
      hit_anything = true;
    }

  return hit_anything;
}

void camera_init(Camera *camera) {
  camera->img_height = (int)((float)camera->img_width / camera->aspect_ratio);

  float viewport_height = 2.0 * tanf(camera->vfov * M_PI / 360.0f) * camera->focal_length;
  float viewport_width = viewport_height * (float)camera->img_width / (float)camera->img_height;

  camera->w = vec3_unit(vec3_sub(camera->look_from, camera->look_to));
  camera->u = vec3_cross(camera->vup, camera->w);
  camera->v = vec3_cross(camera->w, camera->u);

  Vec3 viewport_u = vec3_mul(camera->u, viewport_width);   // scan from left to right
  Vec3 viewport_v = vec3_mul(camera->v, -viewport_height); // scan from top to bottom

  camera->pixel_delta_u = vec3_div(viewport_u, (float)camera->img_width);
  camera->pixel_delta_v = vec3_div(viewport_v, (float)camera->img_height);

  Vec3 viewport_upper_left = vec3_add(camera->look_from, vec3_mul(camera->w, -camera->focal_length),
                                      vec3_mul(viewport_u, -0.5f), vec3_mul(viewport_v, -0.5f));
  camera->pixel00_loc =
      vec3_add(viewport_upper_left, vec3_mul(camera->pixel_delta_u, 0.5f), vec3_mul(camera->pixel_delta_v, 0.5f));

  float dof_radius = camera->focal_length * tanf(camera->dof_angle * M_PI / 360.0f);
  camera->dof_disc_u = vec3_mul(camera->u, dof_radius);
  camera->dof_disc_v = vec3_mul(camera->v, dof_radius);
}

Vec3 camera_ray_color(const Camera *camera, const Ray *ray, const World *world, int depth, PCG32State *rng) {
  if (depth <= 0)
    return (Vec3){0.0f, 0.0f, 0.0f};

  HitRecord hit_record;
  if (hit_objects(world, ray, 1e-3f, INFINITY, &hit_record)) {
    Ray new_ray;
    new_ray.origin = hit_record.p;
    Vec3 scatter_color;

    if (scatter(ray->direction, &hit_record, rng, &new_ray.direction, &scatter_color))
      scatter_color =
          vec3_mul(camera_ray_color(camera, &new_ray, world, depth - 1, rng), scatter_color); // spawn new ray

    return vec3_add(scatter_color, emit(&hit_record));
  }

  // scene background
  return camera->background;

  // old background
  // Vec3 direction = vec3_unit(ray->direction);
  // float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0,1]

  // Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  // Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  // return vec3_lerp(WHITE, BLUE, a);
}

void camera_render(const Camera *camera, const World *world, uint8_t *buffer) {
  for (int j = 0; j < camera->img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", camera->img_height - j);

    int i; // C89 for loop for MSVC OpenMP
#pragma omp parallel for private(i) schedule(static, 1)
    for (i = 0; i < camera->img_width; i++) {
      PCG32State rng;
      pcg32_seed(&rng, 17 + j, 23 + i);

      Vec3 pixel_pos = vec3_add(camera->pixel00_loc, vec3_mul(camera->pixel_delta_u, (float)i),
                                vec3_mul(camera->pixel_delta_v, (float)j));
      Vec3 pixel_color = {0.0f, 0.0f, 0.0f};

      for (int sample = 0; sample < camera->samples_per_pixel; sample++) {
        // square sampling
        // another option: sinc sampling
        // TODO: use 64-bit PRNG to generate 2 numbers at once
        float px = pcg32_f32_between(&rng, -0.5f, 0.5f);
        float py = pcg32_f32_between(&rng, -0.5f, 0.5f);

        Ray ray;
        if (camera->dof_angle > 0.0f) {
          // sample points around camera.look_from like a (thin) lens/aperture
          float a, b;
          for (;;) {
            a = pcg32_f32_between(&rng, -1.0f, 1.0f);
            b = pcg32_f32_between(&rng, -1.0f, 1.0f);
            if (a * a + b * b < 1.0f)
              break;
          }
          ray.origin = vec3_add(camera->look_from, vec3_mul(camera->dof_disc_u, a), vec3_mul(camera->dof_disc_v, b));
        } else {
          ray.origin = camera->look_from;
        }
        ray.direction = vec3_add(pixel_pos, vec3_mul(camera->pixel_delta_u, px), vec3_mul(camera->pixel_delta_v, py),
                                 vec3_neg(ray.origin));

        pixel_color = vec3_add(pixel_color, camera_ray_color(camera, &ray, world, camera->max_depth, &rng));
      }

      pixel_color = vec3_div(pixel_color, (float)camera->samples_per_pixel);
      for (int c = 0; c < 3; c++)
        buffer[(j * camera->img_width + i) * 3 + c] = (int)(256.0f * clamp(sqrtf(pixel_color.x[c]), 0.0f, 0.999f));
    }
  }
  fprintf(stderr, "\nDone\n");
}
