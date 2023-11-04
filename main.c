#include "raytracing.h"
#include "tiff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  World world;
  world.n_materials = 4;
  world.materials = malloc(sizeof(Material) * world.n_materials);
  if (world.materials == NULL) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return 1;
  }
  world.materials[0] = (Material){LAMBERTIAN, {0.8f, 0.8f, 0.0f}};
  world.materials[1] = (Material){DIELECTRIC, {1.0f, 1.0f, 0.8f}, 0.0f, 1.5f};
  world.materials[2] = (Material){METAL, {0.8f, 0.8f, 0.8f}, 0.3f};
  world.materials[3] = (Material){METAL, {0.8f, 0.6f, 0.2f}, 1.0f};

  world.n_spheres = 5;
  world.spheres = malloc(sizeof(Sphere) * world.n_spheres);
  if (world.spheres == NULL) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return 1;
  }
  world.spheres[0] = (Sphere){{0.0f, -100.5f, -1.0f}, 100.0f, world.materials};
  world.spheres[1] = (Sphere){{0.0f, 0.0f, -1.0f}, 0.5f, world.materials + 1};
  world.spheres[2] =
      (Sphere){{0.0f, 0.0f, -1.0f}, -0.4f, world.materials + 1}; // negative radius -> opposite surface normal
  world.spheres[3] = (Sphere){{-1.0f, 0.0f, -1.0f}, 0.5f, world.materials + 2};
  world.spheres[4] = (Sphere){{1.0f, 0.0f, -1.0f}, 0.5f, world.materials + 3};

  Camera camera;
  camera_init(&camera, 16.0f / 9.0f, 400, 10, 5);

  Image8 image = {camera.img_width, camera.img_height, 3, malloc(camera.img_width * camera.img_height * 3)};
  if (image.data == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    return 1;
  }

  camera_render(&camera, &world, image.data);

  FILE *f = fopen("output.tiff", "wb");
  if (f == NULL) {
    fprintf(stderr, "Failed to open file\n");
    return 1;
  }
  write_tiff(f, image);

  return 0;
}
