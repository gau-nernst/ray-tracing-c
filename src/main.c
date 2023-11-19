#include "raytracing.h"
#include "tiff.h"
#include "utils.h"
#include <time.h>

void scene_book1(World *world, Camera *camera) {
  size_t max_spheres = 4 + 22 * 22;
  World_init(world, max_spheres, max_spheres);

  Material mat;

  mat = material(Lambertian_new(texture(Vec3_new(0.5, 0.5, 0.5))));
  MaterialList_append(&world->materials, mat);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, -1000, -1}, 1000, mat)));

  mat = material(Dielectric_new(texture(Vec3_new(1, 1, 1)), 1.5));
  MaterialList_append(&world->materials, mat);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 1, 0}, 1, mat)));

  mat = material(Lambertian_new(texture(Vec3_new(0.4, 0.2, 0.1))));
  MaterialList_append(&world->materials, mat);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){-4, 1, 0}, 1, mat)));

  mat = material(Metal_new(texture(Vec3_new(0.7, 0.6, 0.5)), 0));
  MaterialList_append(&world->materials, mat);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){4, 1, 0}, 1, mat)));

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

        if (choose_material < 0.8f) {
          *color_p = vec3_mul(vec3_rand(&rng), vec3_rand(&rng));
          mat = material(Lambertian_new(texture(color_p)));
        } else if (choose_material < 0.95f) {
          *color_p = vec3_rand_between(&rng, 0.5f, 1);
          mat = material(Metal_new(texture(color_p), pcg32_f32(&rng) * 0.5f));
        } else {
          *color_p = (Vec3){1, 1, 1};
          mat = material(Dielectric_new(texture(color_p), 1.5f));
        }

        MaterialList_append(&world->materials, mat);
        HittableList_append(&world->objects, hittable(Sphere_new(center, radius, mat)));
      }
    }

  BVHNode *bvh = BVHNode_new(world->objects.items, world->objects.size, &rng);
  free(world->objects.items);
  HittableList_init(&world->objects, 1);
  HittableList_append(&world->objects, hittable(bvh));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
  camera->dof_angle = 0.6f;
}

void scene_checker(World *world, Camera *camera) {
  World_init(world, 2, 1);

  Checker *checker_p = Checker_new(0.01f, texture(Vec3_new(0.2, 0.3, 0.1)), texture(Vec3_new(0.9, 0.9, 0.9)));
  Material mat = material(Lambertian_new(texture(checker_p)));
  MaterialList_append(&world->materials, mat);

  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, -10, 0}, 10, mat)));
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 10, 0}, 10, mat)));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
}

void scene_earth(World *world, Camera *camera) {
  World_init(world, 1, 1);

  Material mat = material(Lambertian_new(texture(Image_new("earthmap.jpg"))));
  MaterialList_append(&world->materials, mat);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 0, 0}, 2, mat)));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
}

void scene_perlin(World *world, Camera *camera) {
  World_init(world, 2, 1);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Material mat = material(Lambertian_new(texture(Perlin_new(4.0f, 7, &rng))));
  MaterialList_append(&world->materials, mat);

  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, -1000, 0}, 1000, mat)));
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 2, 0}, 2, mat)));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0.7, 0.8, 1};
  camera->look_from = (Vec3){13, 2, 3};
  camera->look_to = (Vec3){0, 0, 0};
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_simple_light(World *world, Camera *camera) {
  World_init(world, 4, 2);

  PCG32State rng;
  pcg32_seed(&rng, 19, 29);

  Material perlin = material(Lambertian_new(texture(Perlin_new(4.0f, 7, &rng))));
  Material light = material(DiffuseLight_new(texture(Vec3_new(4, 4, 4))));

  MaterialList_append(&world->materials, perlin);
  MaterialList_append(&world->materials, light);

  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, -1000, 0}, 1000, perlin)));
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 2, 0}, 2, perlin)));
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 7, 0}, 2, light)));
  HittableList_append(&world->objects, hittable(Quad_new((Vec3){3, 1, -2}, (Vec3){2, 0, 0}, (Vec3){0, 2, 0}, light)));

  camera->vfov = 20.0f;
  camera->background = (Vec3){0, 0, 0};
  camera->look_from = (Vec3){26, 3, 6};
  camera->look_to = (Vec3){0, 2, 0};
}

void scene_cornell_box(World *world, Camera *camera) {
  World_init(world, 6 + 2, 4);

  Material red = material(Lambertian_new(texture(Vec3_new(0.65, 0.05, 0.05))));
  Material white = material(Lambertian_new(texture(Vec3_new(0.73, 0.73, 0.73))));
  Material green = material(Lambertian_new(texture(Vec3_new(0.12, 0.45, 0.15))));
  Material light = material(DiffuseLight_new(texture(Vec3_new(15, 15, 15))));

  MaterialList_append(&world->materials, red);
  MaterialList_append(&world->materials, white);
  MaterialList_append(&world->materials, green);
  MaterialList_append(&world->materials, light);

  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){555, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, green)));
  HittableList_append(&world->objects, hittable(Quad_new((Vec3){0, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, red)));
  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){343, 554, 332}, (Vec3){-130, 0, 0}, (Vec3){0, 0, -105}, light)));
  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){0, 0, 0}, (Vec3){555, 0, 0}, (Vec3){0, 0, 555}, white)));
  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){555, 555, 555}, (Vec3){-555, 0, 0}, (Vec3){0, 0, -555}, white)));
  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){0, 0, 555}, (Vec3){555, 0, 0}, (Vec3){0, 555, 0}, white)));

  Hittable box1 = hittable(Box_new((Vec3){0, 0, 0}, (Vec3){165, 330, 165}, white));
  box1 = hittable(RotateY_new(box1, 15));
  box1 = hittable(Translate_new(box1, (Vec3){265, 0, 295}));
  HittableList_append(&world->objects, box1);

  Hittable box2 = hittable(Box_new((Vec3){0, 0, 0}, (Vec3){165, 165, 165}, white));
  box2 = hittable(RotateY_new(box2, -18));
  box2 = hittable(Translate_new(box2, (Vec3){130, 0, 65}));
  HittableList_append(&world->objects, box2);

  camera->aspect_ratio = 1.0f;
  camera->background = (Vec3){0, 0, 0};
  camera->vfov = 40.0f;
  camera->look_from = (Vec3){278, 278, -800};
  camera->look_to = (Vec3){278, 278, 0};
}

void scene_book2_final(World *world, Camera *camera) {
  PCG32State rng;
  pcg32_seed(&rng, 19, 29);
  int boxes_per_side = 20;

  World_init(world, boxes_per_side * boxes_per_side + 10, 10);

  Material ground = material(Lambertian_new(texture(Vec3_new(0.48, 0.83, 0.53))));
  MaterialList_append(&world->materials, ground);

  for (int i = 0; i < boxes_per_side; i++)
    for (int j = 0; j < boxes_per_side; j++) {
      float w = 100;
      Vec3 p0 = {-1000 + i * w, 0, -1000 + j * w};
      Vec3 p1 = {-1000 + (i + 1) * w, pcg32_f32_between(&rng, 1, 101), -1000 + (j + 1) * w};
      HittableList_append(&world->objects, hittable(Box_new(p0, p1, ground)));
    }

  Material light = material(DiffuseLight_new(texture(Vec3_new(7, 7, 7))));
  MaterialList_append(&world->materials, light);
  HittableList_append(&world->objects,
                      hittable(Quad_new((Vec3){123, 554, 147}, (Vec3){300, 0, 0}, (Vec3){0, 0, 265}, light)));

  Vec3 center1 = {400, 400, 200};
  Vec3 center2 = {430, 400, 200}; // TODO: moving sphere
  Material sphere_mat = material(Lambertian_new(texture(Vec3_new(0.7, 0.3, 0.1))));
  MaterialList_append(&world->materials, sphere_mat);
  HittableList_append(&world->objects, hittable(Sphere_new(center1, 50, sphere_mat)));

  Material glass = material(Dielectric_new(texture(Vec3_new(1, 1, 1)), 1.5));
  MaterialList_append(&world->materials, glass);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){260, 150, 45}, 50, glass)));

  Material metal = material(Metal_new(texture(Vec3_new(0.8, 0.8, 0.9)), 1.0));
  MaterialList_append(&world->materials, metal);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){0, 150, 145}, 50, metal)));

  // subsurface material
  Hittable boundary = hittable(Sphere_new((Vec3){360, 150, 145}, 70, glass));
  HittableList_append(&world->objects, boundary);
  HittableList_append(&world->objects, hittable(ConstantMedium_new(boundary, 0.2, texture(Vec3_new(0.2, 0.4, 0.9)))));

  // mist
  boundary = hittable(Sphere_new((Vec3){0, 0, 0}, 5000, glass));
  HittableList_append(&world->objects, hittable(ConstantMedium_new(boundary, 0.0001, texture(Vec3_new(1, 1, 1)))));

  Material earth = material(Lambertian_new(texture(Image_new("earthmap.jpg"))));
  MaterialList_append(&world->materials, earth);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){400, 200, 400}, 100, earth)));

  Material perlin = material(Lambertian_new(texture(Perlin_new(0.1, 7, &rng))));
  MaterialList_append(&world->materials, perlin);
  HittableList_append(&world->objects, hittable(Sphere_new((Vec3){220, 280, 300}, 80, perlin)));

  Material white = material(Lambertian_new(texture(Vec3_new(0.73, 0.73, 0.73))));
  MaterialList_append(&world->materials, white);
  int ns = 1000;
  Hittable boxes2 = hittable(HittableList_new(ns));
  for (int i = 0; i < ns; i++)
    HittableList_append(boxes2.ptr, hittable(Sphere_new(vec3_rand_between(&rng, 0, 165), 10, white)));

  boxes2 = hittable(RotateY_new(boxes2, 15));
  boxes2 = hittable(Translate_new(boxes2, (Vec3){-100, 270, 395}));
  HittableList_append(&world->objects, boxes2);

  camera->aspect_ratio = 1.0f;
  camera->background = (Vec3){0, 0, 0};
  camera->vfov = 40.0f;
  camera->look_from = (Vec3){478, 278, -600};
  camera->look_to = (Vec3){278, 278, 0};
}

int main(int argc, char *argv[]) {
  assert(argc > 1);

  World world = {0};
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 400;
  camera.samples_per_pixel = 100;
  camera.max_depth = 10;
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
    fprintf(stderr, "Book 1\n");
    scene_book1(&world, &camera);
    break;
  case 1:
    fprintf(stderr, "Book 2: Checker\n");
    scene_checker(&world, &camera);
    break;
  case 2:
    fprintf(stderr, "Book 2: Earth\n");
    scene_earth(&world, &camera);
    break;
  case 3:
    fprintf(stderr, "Book 2: Perlin noise\n");
    scene_perlin(&world, &camera);
    break;
  case 4:
    fprintf(stderr, "Book 2: Simple light\n");
    scene_simple_light(&world, &camera);
    break;
  case 5:
    fprintf(stderr, "Book 2: Cornell box\n");
    scene_cornell_box(&world, &camera);
    break;
  case 6:
    fprintf(stderr, "Book 2: Final scene\n");
    scene_book2_final(&world, &camera);
    break;
  }
  Camera_init(&camera);

  uint8_t *image = my_malloc(camera.img_width * camera.img_height * 3);

  time_t start, stop;
  time(&start);
  Camera_render(&camera, &world, image);
  time(&stop);
  fprintf(stderr, "Took %ld seconds\n", stop - start);

  FILE *f = fopen("output.tiff", "wb");
  assert((f != NULL) && "Failed to open file");
  write_tiff(f, camera.img_width, camera.img_height, 3, image);

  return 0;
}
