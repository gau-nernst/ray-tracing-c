#include "raytracing.h"
#include "tiff.h"
#include "utils.h"
#include <time.h>

void scene_book1(World *world, Camera *camera) {
  size_t max_spheres = 4 + 22 * 22;
  list_init(&world->spheres, max_spheres);
  list_init(&world->colors, max_spheres);
  list_init(&world->materials, max_spheres);

  list_append(&world->colors, vec3_new(0.5, 0.5, 0.5));
  list_append(&world->materials, lambertian_color_new(list_head(&world->colors)));
  list_append(&world->spheres, sphere_new((Vec3){0, -1000, -1}, 1000, list_head(&world->materials)));

  list_append(&world->colors, vec3_new(1, 1, 1));
  list_append(&world->materials, dielectric_color_new(list_head(&world->colors), 1.5));
  list_append(&world->spheres, sphere_new((Vec3){0, 1, 0}, 1, list_head(&world->materials)));

  list_append(&world->colors, vec3_new(0.4, 0.2, 0.1));
  list_append(&world->materials, lambertian_color_new(list_head(&world->colors)));
  list_append(&world->spheres, sphere_new((Vec3){-4, 1, 0}, 1, list_head(&world->materials)));

  list_append(&world->colors, vec3_new(0.7, 0.6, 0.5));
  list_append(&world->materials, metal_color_new(list_head(&world->colors), 0));
  list_append(&world->spheres, sphere_new((Vec3){4, 1, 0}, 1, list_head(&world->materials)));

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Vec3 ref_point = {4, 0.2, 0};
  float radius = 0.2f;

  for (int a = -11; a < 11; a++)
    for (int b = -11; b < 11; b++) {
      Vec3 center = {
          (float)a + 0.9f * pcg32_f32(&rng),
          radius,
          (float)b + 0.9f * pcg32_f32(&rng),
      };

      if (vec3_length(vec3_sub(center, ref_point)) > 0.9f) {
        float choose_material = pcg32_f32(&rng);

        Vec3 *color_p = my_malloc(sizeof(Vec3));
        Material *material_p;

        if (choose_material < 0.8f) {
          *color_p = vec3_mul(vec3_rand(&rng), vec3_rand(&rng));
          material_p = lambertian_color_new(color_p);
        } else if (choose_material < 0.95f) {
          *color_p = vec3_rand_between(&rng, 0.5f, 1);
          material_p = metal_color_new(color_p, pcg32_f32(&rng) * 0.5f);
        } else {
          *color_p = (Vec3){1, 1, 1};
          material_p = dielectric_color_new(color_p, 1.5f);
        }

        list_append(&world->colors, color_p);
        list_append(&world->materials, material_p);
        list_append(&world->spheres, sphere_new(center, radius, material_p));
      }
    }

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
  camera->dof_angle = 0.6f;
}

// void scene_checker(World *world, Camera *camera) {
//   world->spheres = SphereArray_new(2);
//   world->colors = Vec3Array_new(2);
//   world->checkers = CheckerArray_new(1);
//   world->materials = MaterialArray_new(1);

//   Vec3 *color1 = array_append(&world->colors, vec3(0.2, 0.3, 0.1));
//   Vec3 *color2 = array_append(&world->colors, vec3_full(0.9));

//   Checker *checker_p = array_next(&world->checkers);
//   checker_p->scale = 0.01f;
//   checker_p->even = texture(SOLID, color1);
//   checker_p->odd = texture(SOLID, color2);

//   Material *material_p = array_append(&world->materials, lambertian(texture(CHECKER, checker_p)));
//   array_append(&world->spheres, sphere(vec3(0, -10, 0), 10, material_p));
//   array_append(&world->spheres, sphere(vec3(0, 10, 0), 10, material_p));

//   camera->vfov = 20.0f;
//   camera->background = vec3(0.7, 0.8, 1);
//   camera->look_from = vec3(13, 2, 3);
//   camera->look_to = vec3_zero();
// }

// void scene_earth(World *world, Camera *camera) {
//   world->spheres = SphereArray_new(1);
//   world->materials = MaterialArray_new(1);
//   world->images = ImageArray_new(1);

//   Image *image = array_next(&world->images);
//   image_load(image, "earthmap.jpg");

//   Material *material_p = array_append(&world->materials, lambertian(texture(IMAGE, image)));
//   array_append(&world->spheres, sphere(vec3_zero(), 2, material_p));

//   camera->vfov = 20.0f;
//   camera->background = vec3(0.7, 0.8, 1);
//   camera->look_from = vec3(13, 2, 3);
//   camera->look_to = vec3_zero();
// }

// void scene_perlin(World *world, Camera *camera) {
//   *world = (World){
//       .n_spheres = 2,
//       .n_materials = 1,
//       .n_perlins = 1,
//   };
//   world_malloc(world);

//   PCG32State rng;
//   pcg32_seed(&rng, 19, 29);
//   perlin_init(world->perlins, &rng);
//   world->perlins[0].scale = 4.0f;
//   world->perlins[0].depth = 7;
//   world->materials[0] = (Material){LAMBERTIAN, {PERLIN, .perlin = world->perlins}};

//   world->spheres[0] = sphere(vec3(0, -1000, 0), 1000, world->materials);
//   world->spheres[1] = sphere(vec3(0, 2, 0), 2, world->materials);

//   camera->vfov = 20.0f;
//   camera->background = vec3(0.7, 0.8, 1);
//   camera->look_from = vec3(13, 2, 3);
//   camera->look_to = vec3_zero();
//   camera->dof_angle = 0.0f;
//   camera->focal_length = 10.0f;
// }

// void scene_quads(World *world, Camera *camera) {
//   *world = (World){
//       .n_quads = 5,
//       .n_materials = 5,
//       .n_colors = 5,
//   };
//   world_malloc(world);

//   world->colors[0] = vec3(1, 0.2, 0.2);
//   world->colors[1] = vec3(0.2, 1, 0.2);
//   world->colors[2] = vec3(0.2, 0.2, 1);
//   world->colors[3] = vec3(1, 0.5, 0);
//   world->colors[4] = vec3(0.2, 0.8, 0.8);

//   for (int i = 0; i < world->n_colors; i++)
//     world->materials[i] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + i}};

//   quad_init(world->quads, vec3(-3, -2, 5), vec3(0, 0, -4), vec3(0, 4, 0), world->materials);
//   quad_init(world->quads + 1, vec3(-2, -2, 0), vec3(4, 0, 0), vec3(0, 4, 0), world->materials + 1);
//   quad_init(world->quads + 2, vec3(3, -2, 1), vec3(0, 0, 4), vec3(0, 4, 0), world->materials + 2);
//   quad_init(world->quads + 3, vec3(-2, 3, 1), vec3(4, 0, 0), vec3(0, 0, 4), world->materials + 3);
//   quad_init(world->quads + 4, vec3(-2, -3, 5), vec3(4, 0, 0), vec3(0, 0, -4), world->materials + 4);

//   camera->vfov = 80.0f;
//   camera->background = vec3(0.7, 0.8, 1);
//   camera->look_from = vec3(0, 0, 9);
//   camera->look_to = vec3_zero();
// }

// void scene_simple_light(World *world, Camera *camera) {
//   *world = (World){
//       .n_spheres = 3,
//       .n_quads = 1,
//       .n_materials = 2,
//       .n_colors = 1,
//       .n_perlins = 1,
//   };
//   world_malloc(world);

//   PCG32State rng;
//   pcg32_seed(&rng, 19, 29);
//   perlin_init(world->perlins, &rng);
//   world->perlins[0].scale = 4.0f;
//   world->perlins[0].depth = 7;
//   world->materials[0] = (Material){LAMBERTIAN, {PERLIN, .perlin = world->perlins}};

//   world->colors[0] = vec3_full(4);
//   world->materials[1] = (Material){DIFFUSE_LIGHT, {SOLID, .color = world->colors}};

//   world->spheres[0] = sphere(vec3(0, -1000, 0), 1000, world->materials);
//   world->spheres[1] = sphere(vec3(0, 2, 0), 2, world->materials);
//   world->spheres[2] = sphere(vec3(0, 7, 0), 2, world->materials + 1);
//   quad_init(world->quads, vec3(3, 1, -2), vec3(2, 0, 0), vec3(0, 2, 0), world->materials + 1);

//   camera->vfov = 20.0f;
//   camera->background = vec3_zero();
//   camera->look_from = vec3(26, 3, 6);
//   camera->look_to = vec3(0, 2, 0);
// }

// void scene_cornell_box(World *world, Camera *camera) {
//   *world = (World){
//       .n_quads = 6,
//       .n_colors = 4,
//       .n_materials = 4,
//   };
//   world_malloc(world, 0);

//   world->colors[0] = vec3(0.65, 0.05, 0.05);
//   world->colors[1] = vec3(0.73, 0.73, 0.73);
//   world->colors[2] = vec3(0.12, 0.45, 0.15);
//   for (int i = 0; i < 3; i++)
//     world->materials[i] = (Material){LAMBERTIAN, {SOLID, .color = world->colors + i}};

//   world->colors[3] = vec3_full(15);
//   world->materials[3] = (Material){DIFFUSE_LIGHT, {SOLID, .color = world->colors + 3}};

//   Material *red = world->materials;
//   Material *white = world->materials + 1;
//   Material *green = world->materials + 2;
//   Material *light = world->materials + 3;

//   quad_init(world->quads, vec3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), green);
//   quad_init(world->quads + 1, vec3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), red);
//   quad_init(world->quads + 2, vec3(343, 554, 332), vec3(-130, 0, 0), vec3(0, 0, -105), light);
//   quad_init(world->quads + 3, vec3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), white);
//   quad_init(world->quads + 4, vec3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), white);
//   quad_init(world->quads + 5, vec3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), green);

//   camera->aspect_ratio = 1.0f;
//   camera->background = vec3_zero();
//   camera->vfov = 40.0f;
//   camera->look_from = vec3(278, 278, -800);
//   camera->look_to = vec3(278, 278, 0);
// }

int main(int argc, char *argv[]) {
  assert(argc > 1);

  World world = {0};
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 100;
  camera.max_depth = 50;
  camera.vup = (Vec3){0, 1, 0};
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
    // case 1:
    //   scene_checker(&world, &camera);
    //   break;
    // case 2:
    //   scene_earth(&world, &camera);
    //   break;
    // case 3:
    //   scene_perlin(&world, &camera);
    //   break;
    // case 4:
    //   scene_quads(&world, &camera);
    //   break;
    // case 5:
    //   scene_simple_light(&world, &camera);
    //   break;
    // case 6:
    //   scene_cornell_box(&world, &camera);
    //   break;
  }
  camera_init(&camera);

  uint8_t *image = my_malloc(camera.img_width * camera.img_height * 3);

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
