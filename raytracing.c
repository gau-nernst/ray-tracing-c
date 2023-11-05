#include "raytracing.h"
#define _USE_MATH_DEFINES // for MSVC
#include <math.h>
#include <stdio.h>

Vec3 ray_at(Ray ray, float t) { return vec3vec3_add(ray.origin, vec3float_mul(ray.direction, t)); }

bool hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
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
  hit_record->u = (atan2f(-outward_normal.z, outward_normal.x) + M_PI) * M_1_PI * 0.5;
  hit_record->v = acosf(-outward_normal.y) * M_1_PI;
  hit_record->material = sphere->material;

  return true;
}

bool hit_spheres(const World *world, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  bool hit_anything = false;

  for (int i = 0; i < world->n_spheres; i++)
    if (hit_sphere(world->spheres + i, ray, t_min, t_max, hit_record)) {
      t_max = hit_record->t;
      hit_anything = true;
    }

  return hit_anything;
}

bool scatter_normal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = vec3float_mul(vec3_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
  return false;
}

bool scatter_lambertian(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 new_direction = vec3_add(hit_record->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(new_direction))
    new_direction = hit_record->normal;
  *scattered = new_direction;
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

bool scatter_metal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *scattered = vec3_add(reflect(incident, hit_record->normal),
                        vec3_mul(vec3_rand_unit_vector(rng), hit_record->material->metal_fuzz));
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return vec3_dot(*scattered, hit_record->normal) > 0; // check for degeneration
}

bool scatter_dielectric(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  float eta = hit_record->front_face ? 1.0f / hit_record->material->eta : hit_record->material->eta;
  incident = vec3_unit(incident);

  float cos_theta = fminf(-vec3_dot(incident, hit_record->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_randomf_r(rng)) {
    *scattered = reflect(incident, hit_record->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(incident, vec3_mul(hit_record->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(hit_record->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3_add(r_perp, r_para);
  }
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

bool scatter_default(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = (Vec3){0.0f, 0.0f, 0.0f};
  return false;
}

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  switch (hit_record->material->type) {
  case NORMAL:
    return scatter_normal(incident, hit_record, rng, scattered, color);
  case LAMBERTIAN:
    return scatter_lambertian(incident, hit_record, rng, scattered, color);
  case METAL:
    return scatter_metal(incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return scatter_dielectric(incident, hit_record, rng, scattered, color);
  default:
    return scatter_default(incident, hit_record, rng, scattered, color);
  }
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

Vec3 ray_color(const Ray *ray, const World *world, int depth, PCG32State *rng) {
  if (depth <= 0)
    return (Vec3){0.0f, 0.0f, 0.0f};

  HitRecord hit_record;
  if (hit_spheres(world, ray, 1e-3f, INFINITY, &hit_record)) {
    Ray new_ray;
    new_ray.origin = hit_record.p;
    Vec3 color;

    if (scatter(ray->direction, &hit_record, rng, &new_ray.direction, &color))
      return vec3_mul(ray_color(&new_ray, world, depth - 1, rng), color); // spawn new ray
    return color;
  }

  // scene background
  Vec3 direction = vec3_unit(ray->direction);
  float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0,1]

  Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  return vec3_lerp(WHITE, BLUE, a);
}

void camera_render(Camera *camera, World *world, uint8_t *buffer) {
  for (int j = 0; j < camera->img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", camera->img_height - j);

    int i; // C89 for loop for MSVC OpenMP
#pragma omp parallel for private(i) schedule(static, 1)
    for (i = 0; i < camera->img_width; i++) {
      PCG32State rng;
      pcg32_srandom_r(&rng, 17 + j, 23 + i);

      Vec3 pixel_pos = vec3_add(camera->pixel00_loc, vec3_mul(camera->pixel_delta_u, (float)i),
                                vec3_mul(camera->pixel_delta_v, (float)j));
      Vec3 pixel_color = {0.0f, 0.0f, 0.0f};

      for (int sample = 0; sample < camera->samples_per_pixel; sample++) {
        // square sampling
        // another option: sinc sampling
        // TODO: use 64-bit PRNG to generate 2 numbers at once
        float px = pcg32_randomf_r(&rng) - 0.5f;
        float py = pcg32_randomf_r(&rng) - 0.5f;

        Ray ray;
        if (camera->dof_angle > 0.0f) {
          // sample points around camera.look_from like a (thin) lens/aperture
          float a, b;
          for (;;) {
            a = pcg32_randomf_r(&rng) * 2.0f - 1.0f;
            b = pcg32_randomf_r(&rng) * 2.0f - 1.0f;
            if (a * a + b * b < 1.0f)
              break;
          }
          ray.origin = vec3_add(camera->look_from, vec3_mul(camera->dof_disc_u, a), vec3_mul(camera->dof_disc_v, b));
        } else {
          ray.origin = camera->look_from;
        }
        ray.direction = vec3_add(pixel_pos, vec3_mul(camera->pixel_delta_u, px), vec3_mul(camera->pixel_delta_v, py),
                                 vec3_neg(ray.origin));

        pixel_color = vec3_add(pixel_color, ray_color(&ray, world, camera->max_depth, &rng));
      }

      pixel_color = vec3_div(pixel_color, (float)camera->samples_per_pixel);
      buffer[(j * camera->img_width + i) * 3] = (int)(256.0f * clamp(sqrtf(pixel_color.x), 0.0f, 0.999f));
      buffer[(j * camera->img_width + i) * 3 + 1] = (int)(256.0f * clamp(sqrtf(pixel_color.y), 0.0f, 0.999f));
      buffer[(j * camera->img_width + i) * 3 + 2] = (int)(256.0f * clamp(sqrtf(pixel_color.z), 0.0f, 0.999f));
    }
  }
  fprintf(stderr, "\nDone\n");
}

Vec3 texture_value_solid_color(Texture *texture, float u, float v, Vec3 p) { return texture->color; }
Vec3 texture_value_checker(Texture *texture, float u, float v, Vec3 p) {
  // int x = (int)floorf(p.x / texture->scale);
  // int y = (int)floorf(p.y / texture->scale);
  // int z = (int)floorf(p.z / texture->scale);
  // return texture_value((x + y + z) % 2 ? texture->odd : texture->even, u, v, p);
  int int_u = (int)floorf(u / texture->scale);
  int int_v = (int)floorf(v / texture->scale);
  return texture_value((int_u + int_v) % 2 ? texture->odd : texture->even, u, v, p);
}

Vec3 texture_value(Texture *texture, float u, float v, Vec3 p) {
  switch (texture->type) {
  case SOLID_COLOR:
    return texture_value_solid_color(texture, u, v, p);
  case CHECKER:
    return texture_value_checker(texture, u, v, p);
  default:
    return (Vec3){0.0f, 0.0f, 0.0f};
  }
}
