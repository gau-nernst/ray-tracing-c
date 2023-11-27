#include "raytracing.h"
#include "pcg32.h"
#include <math.h>
#include <stdio.h>

void World_init(World *world, size_t max_objects) {
  HittableList_init(&world->objects, max_objects);
  HittableList_init(&world->lights, max_objects);
}

void Camera_init(Camera *camera) {
  camera->img_height = (int)((float)camera->img_width / camera->aspect_ratio);

  float viewport_height = 2.0f * tanf(camera->vfov * (float)M_PI / 360.0f) * camera->focal_length;
  float viewport_width = viewport_height * (float)camera->img_width / (float)camera->img_height;

  camera->w = vec3_normalize(vec3_sub(camera->look_from, camera->look_to));
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

  float dof_radius = camera->focal_length * tanf(camera->dof_angle * (float)M_PI / 360.0f);
  camera->dof_disc_u = vec3_mul(camera->u, dof_radius);
  camera->dof_disc_v = vec3_mul(camera->v, dof_radius);
}

static Vec3 Camera_ray_color(const Camera *camera, const Ray *ray, const World *world, int depth, PCG32 *rng) {
  if (depth <= 0)
    return VEC3_ZERO;

  HitRecord rec;
  const Hittable *objects = &world->objects.hittable;
  if (objects->vtable->hit(objects, ray, 1e-3f, INFINITY, &rec, rng)) {
    Ray new_ray;
    new_ray.origin = rec.p;
    Vec3 attenuation;
    float sampling_pdf = 0.0f;
    Vec3 emission_color = rec.material->vtable->emit(&rec);

    if (!rec.material->vtable->scatter(&rec, ray->direction, &new_ray.direction, &attenuation, &sampling_pdf, rng))
      return emission_color;

    // mixture pdf
    float prob = camera->lights_sampling_prob;
    if (prob > 0.0f && world->lights.size > 0) {
      const Hittable *lights = &world->lights.hittable;
      if (pcg32_f32(rng) < prob) // change scatter ray towards light source
        new_ray.direction = lights->vtable->rand(lights, rec.p, rng);
      sampling_pdf = (1.0f - prob) * sampling_pdf + prob * lights->vtable->pdf(lights, &new_ray, rng); // update pdf
    }

    // spawn new ray
    Vec3 new_color = Camera_ray_color(camera, &new_ray, world, depth - 1, rng);
    Vec3 scatter_color = vec3_mul(attenuation, new_color);
    if (sampling_pdf > 0.0f) { // TODO: remove this when all materials return sampling_pdf
      float scatter_pdf = rec.material->vtable->scatter_pdf(&rec, ray->direction, new_ray.direction);
      scatter_color = vec3_mul(scatter_color, scatter_pdf / sampling_pdf);
    }
    return vec3_add(scatter_color, emission_color);
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

void Camera_render(const Camera *camera, const World *world, uint8_t *buffer) {
  for (int j = 0; j < camera->img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", camera->img_height - j);

    int i; // C89 for loop for MSVC OpenMP
#pragma omp parallel for private(i) schedule(static, 1)
    for (i = 0; i < camera->img_width; i++) {
      PCG32 rng;
      pcg32_seed(&rng, 17 + j, 23 + i);

      Vec3 pixel_pos = vec3_add(camera->pixel00_loc, vec3_mul(camera->pixel_delta_u, (float)i),
                                vec3_mul(camera->pixel_delta_v, (float)j));
      Vec3 pixel_color = VEC3_ZERO;

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

        pixel_color = vec3_add(pixel_color, Camera_ray_color(camera, &ray, world, camera->max_depth, &rng));
      }

      pixel_color = vec3_div(pixel_color, (float)camera->samples_per_pixel);
      for (int c = 0; c < 3; c++)
        buffer[(j * camera->img_width + i) * 3 + c] = (int)(256.0f * clamp(sqrtf(pixel_color.values[c]), 0.0f, 0.999f));
    }
  }
  fprintf(stderr, "\nDone\n");
}
