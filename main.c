#include "raytracing.h"
#include "tiff.h"
#include <stdio.h>
#include <stdlib.h>

// scene background
Vec3 ray_color(Ray *ray) {
  Vec3 direction = vec3_unit(ray->direction);
  float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0, 1]

  Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  return vec3_add(vec3_mul(WHITE, 1 - a), vec3_mul(BLUE, a));
}

int main(int argc, char *argv[]) {
  float aspect_ratio = 16.0f / 9.0f;
  int img_width = 256;
  int img_height = (int)((float)img_width / aspect_ratio);

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
      vec3_sub(camera_pos, vec3_add(camera_direction, vec3_mul(viewport_u, 0.5f), vec3_mul(viewport_v, 0.5f)));
  Vec3 pixel00_pos = vec3_add(viewport_upper_left, vec3_mul(pixel_delta_u, 0.5f), vec3_mul(pixel_delta_v, 0.5f));

  Image8 image = {img_width, img_height, 3, malloc(img_width * img_height * 3)};
  if (image.data == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    return 1;
  }

  for (int j = 0; j < img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", img_height - j);

    for (int i = 0; i < img_width; i++) {
      Vec3 pixel_pos = vec3_add(pixel00_pos, vec3_mul(pixel_delta_u, (float)i), vec3_mul(pixel_delta_v, (float)j));
      Vec3 ray_direction = vec3_sub(pixel_pos, camera_pos);
      Ray ray = {camera_pos, ray_direction};

      Vec3 pixel_color = ray_color(&ray);

      image.data[(j * img_width + i) * 3] = (int)(255.999 * pixel_color.x);
      image.data[(j * img_width + i) * 3 + 1] = (int)(255.999 * pixel_color.y);
      image.data[(j * img_width + i) * 3 + 2] = (int)(255.999 * pixel_color.z);
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
