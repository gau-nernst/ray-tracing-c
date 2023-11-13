#include "raytracing.h"
#include "tiff.h"
#include "utils.h"
#include <time.h>

void scene_book1(World *world, Camera *camera) {
  size_t max_spheres = 4 + 22 * 22;
  list_init(&world->spheres, max_spheres);
  list_init(&world->materials, max_spheres);

  list_append(&world->materials, Lambertian_new(texture(Vec3_new(0.5, 0.5, 0.5))));
  list_append(&world->spheres, Sphere_new((Vec3){0, -1000, -1}, 1000, list_head(&world->materials)));

  list_append(&world->materials, Dielectric_new(texture(Vec3_new(1, 1, 1)), 1.5));
  list_append(&world->spheres, Sphere_new((Vec3){0, 1, 0}, 1, list_head(&world->materials)));

  list_append(&world->materials, Lambertian_new(texture(Vec3_new(0.4, 0.2, 0.1))));
  list_append(&world->spheres, Sphere_new((Vec3){-4, 1, 0}, 1, list_head(&world->materials)));

  list_append(&world->materials, Metal_new(texture(Vec3_new(0.7, 0.6, 0.5)), 0));
  list_append(&world->spheres, Sphere_new((Vec3){4, 1, 0}, 1, list_head(&world->materials)));

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Vec3 ref_point = {4, 0.2, 0};
  float radius = 0.2f;

  for (int a = -11; a < 11; a++)
    for (int b = -11; b < 11; b++) {
      float choose_material = pcg32_f32(&rng);
      Vec3 center = {
          (float)a + 0.9f * pcg32_f32(&rng),
          radius,
          (float)b + 0.9f * pcg32_f32(&rng),
      };

      if (vec3_length(vec3_sub(center, ref_point)) > 0.9f) {
        Vec3 *color_p = my_malloc(sizeof(Vec3));
        Material *material_p;

        if (choose_material < 0.8f) {
          *color_p = vec3_mul(vec3_rand(&rng), vec3_rand(&rng));
          material_p = Lambertian_new(texture(color_p));
        } else if (choose_material < 0.95f) {
          *color_p = vec3_rand_between(&rng, 0.5f, 1);
          material_p = Metal_new(texture(color_p), pcg32_f32(&rng) * 0.5f);
        } else {
          *color_p = (Vec3){1, 1, 1};
          material_p = Dielectric_new(texture(color_p), 1.5f);
        }

        list_append(&world->materials, material_p);
        list_append(&world->spheres, Sphere_new(center, radius, material_p));
      }
    }

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
  camera->dof_angle = 0.6f;
}

void scene_checker(World *world, Camera *camera) {
  list_init(&world->spheres, 2);
  list_init(&world->materials, 1);

  Checker *checker_p = Checker_new(0.01f, texture(Vec3_new(0.2, 0.3, 0.1)), texture(Vec3_new(0.9, 0.9, 0.9)));
  Material *material_p = Lambertian_new(texture(checker_p));
  list_append(&world->materials, material_p);

  list_append(&world->spheres, Sphere_new((Vec3){0, -10, 0}, 10, material_p));
  list_append(&world->spheres, Sphere_new((Vec3){0, 10, 0}, 10, material_p));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
}

void scene_earth(World *world, Camera *camera) {
  list_init(&world->spheres, 1);
  list_init(&world->materials, 1);

  Material *material_p = Lambertian_new(texture(Image_new("earthmap.jpg")));
  list_append(&world->materials, material_p);
  list_append(&world->spheres, Sphere_new((Vec3){0, 0, 0}, 2, material_p));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
}

void scene_perlin(World *world, Camera *camera) {
  list_init(&world->spheres, 2);
  list_init(&world->materials, 1);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Material *material_p = Lambertian_new(texture(Perlin_new(4.0f, 7, &rng)));
  list_append(&world->materials, material_p);

  list_append(&world->spheres, Sphere_new((Vec3){0, -1000, 0}, 1000, material_p));
  list_append(&world->spheres, Sphere_new((Vec3){0, 2, 0}, 2, material_p));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_simple_light(World *world, Camera *camera) {
  list_init(&world->spheres, 3);
  list_init(&world->quads, 1);
  list_init(&world->materials, 2);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Material *mat_perlin = Lambertian_new(texture(Perlin_new(4.0f, 7, &rng)));
  list_append(&world->materials, mat_perlin);

  Material *mat_light = DiffuseLight_new(texture(Vec3_new(4, 4, 4)));
  list_append(&world->materials, mat_light);

  list_append(&world->spheres, Sphere_new((Vec3){0, -1000, 0}, 1000, mat_perlin));
  list_append(&world->spheres, Sphere_new((Vec3){0, 2, 0}, 2, mat_perlin));
  list_append(&world->spheres, Sphere_new((Vec3){0, 7, 0}, 2, mat_light));
  list_append(&world->quads, Quad_new((Vec3){3, 1, -2}, (Vec3){2, 0, 0}, (Vec3){0, 2, 0}, mat_light));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0, 0, 0};
  camera->look_from = (Vec3){26, 3, 6};
  camera->look_to = (Vec3){0, 2, 0};
}

void scene_cornell_box(World *world, Camera *camera) {
  list_init(&world->quads, 6);
  list_init(&world->materials, 4);

  Material *red = Lambertian_new(texture(Vec3_new(0.65, 0.05, 0.05)));
  Material *white = Lambertian_new(texture(Vec3_new(0.73, 0.73, 0.73)));
  Material *green = Lambertian_new(texture(Vec3_new(0.12, 0.45, 0.15)));
  Material *light = DiffuseLight_new(texture(Vec3_new(15, 15, 15)));

  list_append(&world->materials, red);
  list_append(&world->materials, white);
  list_append(&world->materials, green);
  list_append(&world->materials, light);

  list_append(&world->quads, Quad_new((Vec3){555, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, green));
  list_append(&world->quads, Quad_new((Vec3){0, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, red));
  list_append(&world->quads, Quad_new((Vec3){343, 554, 332}, (Vec3){-130, 0, 0}, (Vec3){0, 0, -105}, light));
  list_append(&world->quads, Quad_new((Vec3){0, 0, 0}, (Vec3){555, 0, 0}, (Vec3){0, 0, 555}, white));
  list_append(&world->quads, Quad_new((Vec3){555, 555, 555}, (Vec3){-555, 0, 0}, (Vec3){0, 0, -555}, white));
  list_append(&world->quads, Quad_new((Vec3){0, 0, 555}, (Vec3){555, 0, 0}, (Vec3){0, 555, 0}, white));

  camera->aspect_ratio = 1.0f;
  camera->background = (Vec3){0, 0, 0};
  camera->vfov = 40.0f;
  camera->look_from = (Vec3){278, 278, -800};
  camera->look_to = (Vec3){278, 278, 0};
}

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
    fprintf(stderr, "Unsupported option. Default to 0\n");
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
    scene_simple_light(&world, &camera);
    break;
  case 5:
    scene_cornell_box(&world, &camera);
    break;
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
