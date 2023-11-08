#include "raytracing.h"
#include "tiff.h"
#include <time.h>

void scene_book1(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 4 + 22 * 22,
      .n_materials = 4 + 22 * 22,
      .n_colors = 4 + 22 * 22,
  };
  world_init(world);

  world->colors[0] = (Vec3){0.5f, 0.5f, 0.5f};
  world->materials[0] = (Material){LAMBERTIAN, {SOLID, .color = world->colors}};
  world->spheres[0] = (Sphere){{0.0f, -1000.0f, -1.0f}, 1000.0f, world->materials};

  world->colors[1] = (Vec3){1.0f, 1.0f, 1.0f};
  world->materials[1] = (Material){DIELECTRIC, {SOLID, .color = world->colors + 1}, 0.0f, 1.5f};
  world->spheres[1] = (Sphere){{0.0f, 1.0f, 0.0f}, 1.0f, world->materials + 1};

  world->colors[2] = (Vec3){0.4f, 0.2f, 0.1f};
  world->materials[2] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + 2}};
  world->spheres[2] = (Sphere){{-4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 2};

  world->colors[3] = (Vec3){0.7f, 0.6f, 0.5f};
  world->materials[3] = (Material){METAL, {SOLID, .color = world->colors + 3}, 0.0f};
  world->spheres[3] = (Sphere){{4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 3};

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Vec3 ref_point = {4.0f, 0.2f, 0.0f};
  size_t index = 4;
  float radius = 0.2f;

  for (int a = -11; a < 11; a++)
    for (int b = -11; b < 11; b++) {
      float choose_material = pcg32_f32(&rng);
      Vec3 center = {(float)a + 0.9f * pcg32_f32(&rng), radius, (float)b + 0.9f * pcg32_f32(&rng)};

      if (vec3_length(vec3_sub(center, ref_point)) > 0.9f) {
        Sphere *sphere = world->spheres + index;
        Material *material = world->materials + index;
        Vec3 *color = world->colors + index;

        *sphere = (Sphere){center, radius, material};
        material->albedo = (Texture){SOLID, .color = color};

        if (choose_material < 0.8f) {
          material->type = LAMBERTIAN;
          *color = vec3_mul(vec3_rand(&rng), vec3_rand(&rng));
        } else if (choose_material < 0.95f) {
          material->type = METAL;
          *color = vec3_rand_between(&rng, 0.5f, 1.0f);
          material->fuzz = pcg32_f32(&rng) * 0.5f;
        } else {
          material->type = DIELECTRIC;
          *color = (Vec3){1.0f, 1.0f, 1.0f};
          material->eta = 1.5f;
        }

        index++;
      }
    }
  world->n_spheres = index;
  world->n_materials = index;
  world->n_colors = index;

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7f, 0.8f, 1.0f};
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->dof_angle = 0.6f;
  camera->focal_length = 10.0f;
}

void scene_checker(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 2,
      .n_materials = 1,
      .n_colors = 2,
      .n_checkers = 1,
  };
  world_init(world);

  world->colors[1] = (Vec3){0.2f, 0.3f, 0.1f};
  world->colors[2] = (Vec3){0.9f, 0.9f, 0.9f};
  world->checkers[0] = (Checker){1e-2f, {SOLID, .color = world->colors}, {SOLID, .color = world->colors + 1}};

  world->materials[0] = (Material){LAMBERTIAN, {CHECKER, .checker = world->checkers}};
  world->spheres[0] = (Sphere){{0.0f, -10.0f, 0.0f}, 10.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 10.0f, 0.0f}, 10.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7f, 0.8f, 1.0f};
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_earth(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 1,
      .n_materials = 1,
      .n_images = 1,
  };
  world_init(world);

  image_load(world->images, "earthmap.jpg");
  world->materials[0] = (Material){LAMBERTIAN, {IMAGE, .image = world->images}};
  world->spheres[0] = (Sphere){{0.0f, 0.0f, 0.0f}, 2.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7f, 0.8f, 1.0f};
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_perlin(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 2,
      .n_materials = 1,
      .n_perlins = 1,
  };
  world_init(world);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);
  perlin_init(world->perlins, &rng);
  world->perlins[0].scale = 4.0f;
  world->perlins[0].depth = 7;
  world->materials[0] = (Material){LAMBERTIAN, {PERLIN, .perlin = world->perlins}};
  world->spheres[0] = (Sphere){{0.0f, -1000.0f, 0.0f}, 1000.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 2.0f, 0.0f}, 2.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7f, 0.8f, 1.0f};
  camera->look_from = (Vec3){13.0f, 2.0f, 3.0f};
  camera->look_to = (Vec3){0.0f, 0.0f, 0.0f};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

int main(int argc, char *argv[]) {
  assert(argc > 1);

  World world;
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 10;
  camera.max_depth = 10;
  camera.vup = (Vec3){0.0f, 1.0f, 0.0f};

  switch (strtol(argv[1], NULL, 10)) {
  default:
  case 0:
    scene_book1(&world, &camera);
    break;
  case 1:
    scene_checker(&world, &camera);
    break;
  case 2:
    scene_earth(&world, &camera);
    break;
  case 3:
    scene_perlin(&world, &camera);
    break;
  }
  camera_init(&camera);

  uint8_t *image;
  try_malloc(image, camera.img_width * camera.img_height * 3);

  time_t start, stop;
  time(&start);
  camera_render(&camera, &world, image);
  time(&stop);
  fprintf(stderr, "Took %ld seconds\n", stop - start);

  FILE *f = fopen("output.tiff", "wb");
  assert((f != NULL) && "Failed to open file");
  write_tiff(f, camera.img_width, camera.img_height, 3, image);

  return 0;
}
