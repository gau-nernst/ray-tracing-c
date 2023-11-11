#include "raytracing.h"
#include "tiff.h"
#include <time.h>

void scene_book1(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 4 + 22 * 22,
      .n_materials = 4 + 22 * 22,
      .n_colors = 4 + 22 * 22,
  };
  world_malloc(world);

  world->colors[0] = vec3_full(0.5f);
  world->materials[0] = (Material){LAMBERTIAN, {SOLID, .color = world->colors}};
  world->spheres[0] = (Sphere){{0.0f, -1000.0f, -1.0f}, 1000.0f, world->materials};

  world->colors[1] = vec3_full(1.0f);
  world->materials[1] = (Material){DIELECTRIC, {SOLID, .color = world->colors + 1}, 0.0f, 1.5f};
  world->spheres[1] = (Sphere){{0.0f, 1.0f, 0.0f}, 1.0f, world->materials + 1};

  world->colors[2] = vec3(0.4f, 0.2f, 0.1f);
  world->materials[2] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + 2}};
  world->spheres[2] = (Sphere){{-4.0f, 1.0f, 0.0f}, 1.0f, world->materials + 2};

  world->colors[3] = vec3(0.7f, 0.6f, 0.5f);
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
          *color = vec3_full(1.0f);
          material->eta = 1.5f;
        }

        index++;
      }
    }
  world->n_spheres = index;
  world->n_materials = index;
  world->n_colors = index;

  camera->vfov = 20.0f;
  camera->background = vec3(0.7f, 0.8f, 1.0f);
  camera->look_from = vec3(13.0f, 2.0f, 3.0f);
  camera->look_to = vec3_zero();
  camera->dof_angle = 0.6f;
}

void scene_checker(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 2,
      .n_materials = 1,
      .n_colors = 2,
      .n_checkers = 1,
  };
  world_malloc(world);

  world->colors[0] = vec3(0.2f, 0.3f, 0.1f);
  world->colors[1] = vec3(0.9f, 0.9f, 0.9f);
  world->checkers[0] = (Checker){1e-2f, {SOLID, .color = world->colors}, {SOLID, .color = world->colors + 1}};

  world->materials[0] = (Material){LAMBERTIAN, {CHECKER, .checker = world->checkers}};
  world->spheres[0] = (Sphere){{0.0f, -10.0f, 0.0f}, 10.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 10.0f, 0.0f}, 10.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = vec3(0.7f, 0.8f, 1.0f);
  camera->look_from = vec3(13.0f, 2.0f, 3.0f);
  camera->look_to = vec3_zero();
}

void scene_earth(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 1,
      .n_materials = 1,
      .n_images = 1,
  };
  world_malloc(world);

  image_load(world->images, "earthmap.jpg");
  world->materials[0] = (Material){LAMBERTIAN, {IMAGE, .image = world->images}};
  world->spheres[0] = (Sphere){vec3_zero(), 2.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = vec3(0.7f, 0.8f, 1.0f);
  camera->look_from = vec3(13.0f, 2.0f, 3.0f);
  camera->look_to = vec3_zero();
}

void scene_perlin(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 2,
      .n_materials = 1,
      .n_perlins = 1,
  };
  world_malloc(world);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);
  perlin_init(world->perlins, &rng);
  world->perlins[0].scale = 4.0f;
  world->perlins[0].depth = 7;
  world->materials[0] = (Material){LAMBERTIAN, {PERLIN, .perlin = world->perlins}};

  world->spheres[0] = (Sphere){{0.0f, -1000.0f, 0.0f}, 1000.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 2.0f, 0.0f}, 2.0f, world->materials};

  camera->vfov = 20.0f;
  camera->background = vec3(0.7f, 0.8f, 1.0f);
  camera->look_from = vec3(13.0f, 2.0f, 3.0f);
  camera->look_to = vec3_zero();
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_quads(World *world, Camera *camera) {
  *world = (World){
      .n_quads = 5,
      .n_materials = 5,
      .n_colors = 5,
  };
  world_malloc(world);

  world->colors[0] = vec3(1.0f, 0.2f, 0.2f);
  world->colors[1] = vec3(0.2f, 1.0f, 0.2f);
  world->colors[2] = vec3(0.2f, 0.2f, 1.0f);
  world->colors[3] = vec3(1.0f, 0.5f, 0.0f);
  world->colors[4] = vec3(0.2f, 0.8f, 0.8f);

  world->quads[0] = (Quad){{-3.0f, -2.0f, 5.0f}, {0.0f, 0.0f, -4.0f}, {0.0f, 4.0f, 0.0f}};
  world->quads[1] = (Quad){{-2.0f, -2.0f, 0.0f}, {4.0f, 0.0f, 0.0f}, {0.0f, 4.0f, 0.0f}};
  world->quads[2] = (Quad){{3.0f, -2.0f, 1.0f}, {0.0f, 0.0f, 4.0f}, {0.0f, 4.0f, 0.0f}};
  world->quads[3] = (Quad){{-2.0f, 3.0f, 1.0f}, {4.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 4.0f}};
  world->quads[4] = (Quad){{-2.0f, -3.0f, 5.0f}, {4.0f, 0.0f, 0.0f}, {0.0f, 0.0f, -4.0f}};

  for (int i = 0; i < world->n_colors; i++) {
    world->materials[i] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + i}};
    world->quads[i].material = world->materials + i;
    quad_init(world->quads + i);
  }

  camera->vfov = 80.0f;
  camera->background = vec3(0.7f, 0.8f, 1.0f);
  camera->look_from = vec3(0.0f, 0.0f, 9.0f);
  camera->look_to = vec3_zero();
}

void scene_simple_light(World *world, Camera *camera) {
  *world = (World){
      .n_spheres = 3,
      .n_quads = 1,
      .n_materials = 2,
      .n_colors = 1,
      .n_perlins = 1,
  };
  world_malloc(world);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);
  perlin_init(world->perlins, &rng);
  world->perlins[0].scale = 4.0f;
  world->perlins[0].depth = 7;
  world->materials[0] = (Material){LAMBERTIAN, {PERLIN, .perlin = world->perlins}};

  world->colors[0] = vec3(4.0f, 4.0f, 4.0f);
  world->materials[1] = (Material){DIFFUSE_LIGHT, {SOLID, .color = world->colors}};

  world->spheres[0] = (Sphere){{0.0f, -1000.0f, 0.0f}, 1000.0f, world->materials};
  world->spheres[1] = (Sphere){{0.0f, 2.0f, 0.0f}, 2.0f, world->materials};
  world->spheres[2] = (Sphere){{0.0f, 7.0f, 0.0f}, 2.0f, world->materials + 1};
  world->quads[0] =
      (Quad){{3.0f, 1.0f, -2.0f}, {2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}, .material = world->materials + 1};
  quad_init(world->quads);

  camera->vfov = 20.0f;
  camera->background = vec3_zero();
  camera->look_from = vec3(26.0f, 3.0f, 6.0f);
  camera->look_to = vec3(0.0f, 2.0f, 0.0f);
}

void scene_cornell_box(World *world, Camera *camera) {
  *world = (World){
      .n_quads = 6,
      .n_colors = 4,
      .n_materials = 4,
  };
  world_malloc(world);

  world->colors[0] = vec3(0.65f, 0.05f, 0.05f);
  world->colors[1] = vec3(0.73f, 0.73f, 0.73f);
  world->colors[2] = vec3(0.12f, 0.45f, 0.15f);
  world->colors[3] = vec3(15.0f, 15.0f, 15.0f);

  for (int i = 0; i < 3; i++)
    world->materials[i] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + i}};
  world->materials[3] = (Material){DIFFUSE_LIGHT, {SOLID, .color = world->colors + 3}};

  Material *red = world->materials;
  Material *white = world->materials + 1;
  Material *green = world->materials + 2;
  Material *light = world->materials + 3;

  world->quads[0] = (Quad){{555.0f, 0.0f, 0.0f}, {0.0f, 555.0f, 0.0f}, {0.0f, 0.0f, 555.0f}, .material = green};
  world->quads[1] = (Quad){{0.0f, 0.0f, 0.0f}, {0.0f, 555.0f, 0.0f}, {0.0f, 0.0f, 555.0f}, .material = red};
  world->quads[2] = (Quad){{343.0f, 554.0f, 332.0f}, {-130.0f, 0.0f, 0.0f}, {0.0f, 0.0f, -105.0f}, .material = light};
  world->quads[3] = (Quad){{0.0f, 0.0f, 0.0f}, {555.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 555.0f}, .material = white};
  world->quads[4] = (Quad){{555.0f, 555.0f, 555.0f}, {-555.0f, 0.0f, 0.0f}, {0.0f, 0.0f, -555.0f}, .material = white};
  world->quads[5] = (Quad){{0.0f, 0.0f, 555.0f}, {555.0f, 0.0f, 0.0f}, {0.0f, 555.0f, 0.0f}, .material = white};
  for (int i = 0; i < 6; i++)
    quad_init(world->quads + i);

  camera->aspect_ratio = 1.0f;
  camera->background = vec3_zero();
  camera->vfov = 40.0f;
  camera->look_from = vec3(278.0f, 278.0f, -800.0f);
  camera->look_to = vec3(278.0f, 278.0f, 0.0f);
}

int main(int argc, char *argv[]) {
  assert(argc > 1);

  World world;
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 100;
  camera.max_depth = 50;
  camera.vup = vec3(0.0f, 1.0f, 0.0f);
  camera.dof_angle = 0.0f;
  camera.focal_length = 10.0f;

  if (argc > 3)
    camera.img_width = strtol(argv[2], NULL, 10);
  if (argc > 4)
    camera.samples_per_pixel = strtol(argv[3], NULL, 10);

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
  case 4:
    scene_quads(&world, &camera);
    break;
  case 5:
    scene_simple_light(&world, &camera);
  case 6:
    scene_cornell_box(&world, &camera);
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
