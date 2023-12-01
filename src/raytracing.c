#include "raytracing.h"
#include "material.h"
#include "pcg32.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

void World_init(World *world, size_t max_objects) {
  HittableList_init(&world->objects, max_objects);
  HittableList_init(&world->lights, max_objects);
}

void Camera_init(Camera *self, Vec3 background, float vfov, float focal_length, Vec3 look_from, Vec3 look_to,
                 float dof_angle) {
  self->background = background;
  self->vfov = vfov;
  self->focal_length = focal_length;

  self->look_from = look_from;
  self->look_to = look_to;
  self->vup = vec3(0, 1, 0);
  self->w = vec3_normalize(vec3_sub(look_from, look_to));
  self->u = vec3_cross(self->vup, self->w);
  self->v = vec3_cross(self->w, self->u);

  self->dof_angle = dof_angle;
  float dof_radius = focal_length * tanf(dof_angle * (float)M_PI / 360.0f);
  self->dof_disc_u = vec3_mul(self->u, dof_radius);
  self->dof_disc_v = vec3_mul(self->v, dof_radius);
}

void Renderer_init(Renderer *self, int img_width, float aspect_ratio, int samples_per_pixel, int max_depth,
                   float lights_sampling_prob) {
  self->img_width = img_width;
  self->img_height = (int)((float)img_width / aspect_ratio);
  self->samples_per_pixel = samples_per_pixel;
  self->max_depth = max_depth;
  self->lights_sampling_prob = lights_sampling_prob;
}

static Vec3 ray_color(const Ray *ray, const World *world, int depth, float lights_sampling_prob, Vec3 background,
                      PCG32 *rng) {
  if (depth <= 0)
    return VEC3_ZERO;

  HitRecord rec;
  const Hittable *objects = &world->objects.hittable;
  if (objects->vtable->hit(objects, ray, 1e-3f, INFINITY, &rec, rng)) {
    Ray r_out = {rec.p};
    Vec3 albedo;
    bool skip_pdf;
    Vec3 emission_color = Material_emit(&rec);

    if (!Material_scatter(&rec, ray->direction, &r_out.direction, &albedo, &skip_pdf, rng))
      return emission_color;

    if (skip_pdf || lights_sampling_prob == 0.0f || world->lights.size == 0) {
      Vec3 scatter_color = vec3_mul(albedo, ray_color(&r_out, world, depth - 1, lights_sampling_prob, background, rng));
      return vec3_add(emission_color, scatter_color);
    }

    // mixture pdf: change scatter ray towards light source
    const Hittable *lights = &world->lights.hittable;
    if (pcg32_f32(rng) < lights_sampling_prob)
      r_out.direction = lights->vtable->rand(lights, rec.p, rng);

    float scatter_pdf = Material_scatter_pdf(rec.material, rec.normal, ray->direction, r_out.direction);
    float sampling_pdf =
        (1.0f - lights_sampling_prob) * scatter_pdf + lights_sampling_prob * lights->vtable->pdf(lights, &r_out, rng);

    Vec3 scatter_color = vec3_mul(albedo, ray_color(&r_out, world, depth - 1, lights_sampling_prob, background, rng),
                                  scatter_pdf / sampling_pdf);
    return vec3_add(emission_color, scatter_color);
  }

  // scene background
  return background;

  // old background
  // Vec3 direction = vec3_unit(ray->direction);
  // float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0,1]

  // Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  // Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  // return vec3_lerp(WHITE, BLUE, a);
}

void render(const Renderer *renderer, const Camera *camera, const World *world, uint8_t *buffer) {
  float viewport_height = 2.0f * tanf(camera->vfov * (float)M_PI / 360.0f) * camera->focal_length;
  float viewport_width = viewport_height * (float)renderer->img_width / (float)renderer->img_height;

  Vec3 viewport_u = vec3_mul(camera->u, viewport_width);   // scan from left to right
  Vec3 viewport_v = vec3_mul(camera->v, -viewport_height); // scan from top to bottom
  Vec3 viewport_upper_left = vec3_add(camera->look_from, vec3_mul(camera->w, -camera->focal_length),
                                      vec3_mul(viewport_u, -0.5f), vec3_mul(viewport_v, -0.5f));

  Vec3 pixel_delta_u = vec3_div(viewport_u, (float)renderer->img_width);
  Vec3 pixel_delta_v = vec3_div(viewport_v, (float)renderer->img_height);
  Vec3 pixel00_loc = vec3_add(viewport_upper_left, vec3_mul(pixel_delta_u, 0.5f), vec3_mul(pixel_delta_v, 0.5f));

  for (int j = 0; j < renderer->img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", renderer->img_height - j);

    int i; // C89 for loop for MSVC OpenMP
#pragma omp parallel for private(i) schedule(static, 1)
    for (i = 0; i < renderer->img_width; i++) {
      PCG32 rng;
      pcg32_seed(&rng, 17 + j, 23 + i);

      Vec3 pixel_pos = vec3_add(pixel00_loc, vec3_mul(pixel_delta_u, (float)i), vec3_mul(pixel_delta_v, (float)j));
      Vec3 pixel_color = VEC3_ZERO;

      for (int sample = 0; sample < renderer->samples_per_pixel; sample++) {
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
        ray.direction =
            vec3_add(pixel_pos, vec3_mul(pixel_delta_u, px), vec3_mul(pixel_delta_v, py), vec3_neg(ray.origin));

        pixel_color = vec3_add(pixel_color, ray_color(&ray, world, renderer->max_depth, renderer->lights_sampling_prob,
                                                      camera->background, &rng));
      }

      for (int c = 0; c < 3; c++) {
        float value = pixel_color.values[c];
        value = clamp(sqrtf(value / renderer->samples_per_pixel), 0.0f, 0.999f);
        buffer[(j * renderer->img_width + i) * 3 + c] = (int)(256.0f * value);
      }
    }
  }
  fprintf(stderr, "\nDone\n");
}
