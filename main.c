#include "raytracing.h"
#include "tiff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int scene1(World *world, Camera *camera) {
  world->n_spheres = 4 + 22 * 22;
  world->n_materials = 4 + 22 * 22;
  world->spheres = malloc(sizeof(Sphere) * world->n_spheres);
  world->materials = malloc(sizeof(Material) * world->n_materials);
  if (world->spheres == NULL || world->materials == NULL) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return 1;
  }

  world->materials[0] = (Material){LAMBERTIAN, {SOLID_COLOR, {0.5f, 0.5f, 0.5f}}};
  world->materials[1] = (Material){DIELECTRIC, {SOLID_COLOR, {1.0f, 1.0f, 1.0f}}, 0.0f, 1.5f};
  world->materials[2] = (Material){LAMBERTIAN, {SOLID_COLOR, {0.4f, 0.2f, 0.1f}}};
  world->materials[3] = (Material){METAL, {SOLID_COLOR, {0.7f, 0.6f, 0.5f}}, 0.0f};

  world->spheres[0] = (Sphere){{0.0f, -1000.0f, -1.0f}, 1000.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 1.0f, 0.0f}, 1.0f, world->materials + 1};
  world->spheres[2] = (Sphere){{-4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 2};
  world->spheres[3] = (Sphere){{4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 3};

  PCG32State rng;
  pcg32_srandom_r(&rng, 19, 29);

  Vec3 ref_point = {4.0f, 0.2f, 0.0f};
  size_t index = 4;
  float radius = 0.2f;

  for (int a = -11; a < 11; a++)
    for (int b = -11; b < 11; b++) {
      float choose_material = pcg32_randomf_r(&rng);
      Vec3 center = {(float)a + 0.9f * pcg32_randomf_r(&rng), radius, (float)b + 0.9f * pcg32_randomf_r(&rng)};

      if (vec3_length(vec3_sub(center, ref_point)) > 0.9f) {
        Material *material = world->materials + index;
        if (choose_material < 0.8f) {
          material->type = LAMBERTIAN;
          material->albedo = (Texture){SOLID_COLOR, vec3_mul(vec3_rand(&rng), vec3_rand(&rng))};
        } else if (choose_material < 0.95f) {
          material->type = METAL;
          material->albedo = (Texture){SOLID_COLOR, vec3_rand_between(0.5f, 1.0f, &rng)};
          material->metal_fuzz = pcg32_randomf_r(&rng) / 2.0f;
        } else {
          material->type = DIELECTRIC;
          material->albedo = (Texture){SOLID_COLOR, {1.0f, 1.0f, 1.0f}};
          material->eta = 1.5f;
        }

        world->spheres[index].center = center;
        world->spheres[index].radius = radius;
        world->spheres[index].material = material;
        index++;
      }
    }
  world->n_spheres = index;
  world->n_materials = index;

  camera->vfov = 20.0f;
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->vup = (Vec3){0.0f, 1.0f, 0.0f};
  camera->dof_angle = 0.6f;
  camera->focal_length = 10.0f;

  return 0;
}

int scene2(World *world, Camera *camera) {
  Texture *textures = malloc(sizeof(Texture) * 2);
  world->materials = malloc(sizeof(Material) * world->n_materials);
  world->spheres = malloc(sizeof(Sphere) * world->n_spheres);
  if (textures == NULL || world->spheres == NULL || world->materials == NULL) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return 1;
  }

  textures[0] = (Texture){SOLID_COLOR, {0.2f, 0.3f, 0.1f}};
  textures[1] = (Texture){SOLID_COLOR, {0.9f, 0.9f, 0.9f}};

  world->n_materials = 1;
  world->materials[0] = (Material){LAMBERTIAN, {CHECKER, {0.0f, 0.0f, 0.0f}, 1e-2f, textures, textures + 1}};

  world->n_spheres = 2;
  world->spheres[0] = (Sphere){{0.0f, -10.0f, 0.0f}, 10.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 10.0f, 0.0f}, 10.0f, world->materials};

  camera->vfov = 20.0f;
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->vup = (Vec3){0.0f, 1.0f, 0.0f};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;

  return 0;
}

int main(int argc, char *argv[]) {
  World world;
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 10;
  camera.max_depth = 10;

  if (scene2(&world, &camera))
    return 1;
  camera_init(&camera);

  Image8 image = {camera.img_width, camera.img_height, 3, malloc(camera.img_width * camera.img_height * 3)};
  if (image.data == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    return 1;
  }

  time_t start, stop;
  time(&start);
  camera_render(&camera, &world, image.data);
  time(&stop);
  fprintf(stderr, "Took %ld seconds\n", stop - start);

  FILE *f = fopen("output.tiff", "wb");
  if (f == NULL) {
    fprintf(stderr, "Failed to open file\n");
    return 1;
  }
  write_tiff(f, image);

  return 0;
}
