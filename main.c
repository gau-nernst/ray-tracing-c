#include "raytracing.h"
#include "tiff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#define clamp(x, lo, hi) min(max(x, lo), hi)

// scene background
Vec3 ray_color(const Ray *ray, const Sphere *spheres, int n) {
  HitRecord hit_record;
  if (hit_spheres(spheres, n, ray, 0.0f, 100000.0f, &hit_record)) {
    return vec3_mul(vec3_add(hit_record.normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
  }

  Vec3 direction = vec3_unit(ray->direction);
  float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0,1]

  Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  return vec3_add(vec3_mul(WHITE, 1 - a), vec3_mul(BLUE, a));
}

int main(int argc, char *argv[]) {
  float aspect_ratio = 16.0f / 9.0f;
  int img_width = 400;
  int img_height = (int)((float)img_width / aspect_ratio);

  Sphere *spheres = malloc(sizeof(Sphere) * 2);
  if (spheres == NULL) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return 1;
  }
  spheres[0] = (Sphere){{0.0f, 0.0f, -1.0f}, 0.5f};
  spheres[1] = (Sphere){{0.0f, -100.5f, -1.0f}, 100.0f};

  int samples_per_pixel = 10;
  float focal_length = 1.0f;
  float viewport_height = 2.0f;
  float viewport_width = viewport_height * (float)img_width / (float)img_height;
  Vec3 camera_pos = {0.0f, 0.0f, 0.0f};

  Vec3 viewport_u = {viewport_width, 0.0f, 0.0f};   // scan from left to right
  Vec3 viewport_v = {0.0f, -viewport_height, 0.0f}; // scan from top to bottom

  Vec3 pixel_delta_u = vec3_mul(viewport_u, 1.0f / (float)img_width);
  Vec3 pixel_delta_v = vec3_mul(viewport_v, 1.0f / (float)img_height);

  Vec3 camera_direction = {0.0f, 0.0f, focal_length};
  Vec3 viewport_upper_left =
      vec3_add(camera_pos, vec3_neg(camera_direction), vec3_mul(viewport_u, -0.5f), vec3_mul(viewport_v, -0.5f));
  Vec3 pixel00_pos = vec3_add(viewport_upper_left, vec3_mul(pixel_delta_u, 0.5f), vec3_mul(pixel_delta_v, 0.5f));

  Image8 image = {img_width, img_height, 3, malloc(img_width * img_height * 3)};
  if (image.data == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    return 1;
  }

  PCG32State pcg32_state;
  pcg32_srandom_r(&pcg32_state, 10, 20);

  for (int j = 0; j < img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", img_height - j);

    for (int i = 0; i < img_width; i++) {
      Vec3 pixel_pos = vec3_add(pixel00_pos, vec3_mul(pixel_delta_u, (float)i), vec3_mul(pixel_delta_v, (float)j));
      Vec3 pixel_color = {0.0f, 0.0f, 0.0f};

      for (int sample = 0; sample < samples_per_pixel; sample++) {
        // square sample
        // TODO: use 64-bit PRNG to generate 2 numbers at once
        float px = pcg32_randomf_r(&pcg32_state) - 0.5f;
        float py = pcg32_randomf_r(&pcg32_state) - 0.5f;
        Vec3 ray_direction =
            vec3_add(pixel_pos, vec3_mul(pixel_delta_u, px), vec3_mul(pixel_delta_v, py), vec3_neg(camera_pos));
        Ray ray = {camera_pos, ray_direction};
        pixel_color = vec3_add(pixel_color, ray_color(&ray, spheres, 2));
      }

      pixel_color = vec3_mul(pixel_color, 1.0f / (float)samples_per_pixel);
      image.data[(j * img_width + i) * 3] = (int)(256.0f * clamp(pixel_color.x, 0.0f, 0.999f));
      image.data[(j * img_width + i) * 3 + 1] = (int)(256.0f * clamp(pixel_color.y, 0.0f, 0.999f));
      image.data[(j * img_width + i) * 3 + 2] = (int)(256.0f * clamp(pixel_color.z, 0.0f, 0.999f));
    }
  }
  fprintf(stderr, "\nDone\n");

  FILE *f = fopen("output.tiff", "wb");
  if (f == NULL) {
    fprintf(stderr, "Failed to open file\n");
    return 1;
  }
  write_tiff(f, image);

  return 0;
}
