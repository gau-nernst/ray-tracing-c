#include "raytracing.h"
#include "tiff.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define try_malloc(ptr, sz) assert(((ptr = malloc(sz)) != NULL) && "Failed to allocate memory")

void world_init(World *world) {
  try_malloc(world->spheres, sizeof(Sphere) * world->n_spheres);
  try_malloc(world->materials, sizeof(Material) * world->n_materials);
  try_malloc(world->textures, sizeof(Texture) * world->n_textures);
}

void scene1(World *world, Camera *camera) {
  world->n_spheres = 4 + 22 * 22;
  world->n_materials = 4 + 22 * 22;
  world->n_textures = 4 + 22 * 22;
  world_init(world);

  world->textures[0] = (Texture){SOLID, .color = {0.5f, 0.5f, 0.5f}};
  world->materials[0] = (Material){LAMBERTIAN, world->textures};
  world->spheres[0] = (Sphere){{0.0f, -1000.0f, -1.0f}, 1000.0f, world->materials};

  world->textures[1] = (Texture){SOLID, .color = {1.0f, 1.0f, 1.0f}};
  world->materials[1] = (Material){DIELECTRIC, world->textures + 1, 0.0f, 1.5f};
  world->spheres[1] = (Sphere){{0.0f, 1.0f, 0.0f}, 1.0f, world->materials + 1};

  world->textures[2] = (Texture){SOLID, .color = {0.4f, 0.2f, 0.1f}};
  world->materials[2] = (Material){LAMBERTIAN, world->textures + 2};
  world->spheres[2] = (Sphere){{-4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 2};

  world->textures[3] = (Texture){SOLID, .color = {0.7f, 0.6f, 0.5f}};
  world->materials[3] = (Material){METAL, world->textures + 3, 0.0f};
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
        Sphere *sphere = world->spheres + index;
        Material *material = world->materials + index;
        Texture *texture = world->textures + index;

        *sphere = (Sphere){center, radius, material};
        material->albedo = texture;

        if (choose_material < 0.8f) {
          material->type = LAMBERTIAN;
          *texture = (Texture){SOLID, .color = vec3_mul(vec3_rand(&rng), vec3_rand(&rng))};
        } else if (choose_material < 0.95f) {
          material->type = METAL;
          *texture = (Texture){SOLID, .color = vec3_rand_between(0.5f, 1.0f, &rng)};
          material->fuzz = pcg32_randomf_r(&rng) / 2.0f;
        } else {
          material->type = DIELECTRIC;
          *texture = (Texture){SOLID, .color = {1.0f, 1.0f, 1.0f}};
          material->eta = 1.5f;
        }

        index++;
      }
    }
  world->n_spheres = index;
  world->n_materials = index;
  world->n_textures = index;

  camera->vfov = 20.0f;
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->vup = (Vec3){0.0f, 1.0f, 0.0f};
  camera->dof_angle = 0.6f;
  camera->focal_length = 10.0f;
}

void scene2(World *world, Camera *camera) {
  world->n_textures = 3;
  world->n_materials = 1;
  world->n_spheres = 2;
  world_init(world);

  world->textures[0] = (Texture){SOLID, .color = {0.2f, 0.3f, 0.1f}};
  world->textures[1] = (Texture){SOLID, .color = {0.9f, 0.9f, 0.9f}};
  world->textures[2] = (Texture){CHECKER, .scale = 1e-2f, .even = world->textures, .odd = world->textures + 1};

  world->materials[0] = (Material){LAMBERTIAN, world->textures + 2};
  world->spheres[0] = (Sphere){{0.0f, -10.0f, 0.0f}, 10.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 10.0f, 0.0f}, 10.0f, world->materials};

  camera->vfov = 20.0f;
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->vup = (Vec3){0.0f, 1.0f, 0.0f};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

int main(int argc, char *argv[]) {
  World world;
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 10;
  camera.max_depth = 10;

  scene2(&world, &camera);
  camera_init(&camera);

  Image8 image = {camera.img_width, camera.img_height, 3};
  try_malloc(image.data, camera.img_width * camera.img_height * 3);

  time_t start, stop;
  time(&start);
  camera_render(&camera, &world, image.data);
  time(&stop);
  fprintf(stderr, "Took %ld seconds\n", stop - start);

  FILE *f = fopen("output.tiff", "wb");
  assert((f != NULL) && "Failed to open file");
  write_tiff(f, image);

  return 0;
}
