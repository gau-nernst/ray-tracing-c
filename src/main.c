#include "hittable.h"
#include "material.h"
#include "raytracing.h"
#include "texture.h"
#include "tiff.h"
#include "utils.h"
#include <time.h>

void scene_metal_and_lambertian(World *world, Camera *camera) {
  World_init(world, 4);

  Material *mat;

  mat = Lambertian_new(Solid_new(vec3(0.8, 0.8, 0.0)));
  HittableList_append(&world->objects, Sphere_new(vec3(0, -100.5, -1), 100, mat));

  mat = Lambertian_new(Solid_new(vec3(0.7, 0.3, 0.3)));
  HittableList_append(&world->objects, Sphere_new(vec3(0, 0, -1), 0.5, mat));

  mat = Metal_new(Solid_new(vec3(0.8, 0.8, 0.8)), 0.3);
  HittableList_append(&world->objects, Sphere_new(vec3(-1, 0, -1), 0.5, mat));

  mat = Metal_new(Solid_new(vec3(0.8, 0.6, 0.2)), 1.0);
  HittableList_append(&world->objects, Sphere_new(vec3(1, 0, -1), 0.5, mat));

  camera->vfov = 90.0f;
  camera->background = vec3(0.7, 0.8, 1);
  camera->look_from = VEC3_ZERO;
  camera->look_to = vec3(0, 0, -1);
}

void scene_book1_final(World *world, Camera *camera) {
  World_init(world, 4 + 22 * 22);

  Material *mat;

  mat = Lambertian_new(Solid_new(vec3(0.5, 0.5, 0.5)));
  HittableList_append(&world->objects, Sphere_new(vec3(0, -1000, -1), 1000, mat));

  mat = Dielectric_new(1.5);
  HittableList_append(&world->objects, Sphere_new(vec3(0, 1, 0), 1, mat));

  mat = Lambertian_new(Solid_new(vec3(0.4, 0.2, 0.1)));
  HittableList_append(&world->objects, Sphere_new(vec3(-4, 1, 0), 1, mat));

  mat = Metal_new(Solid_new(vec3(0.7, 0.6, 0.5)), 0);
  HittableList_append(&world->objects, Sphere_new(vec3(4, 1, 0), 1, mat));

  PCG32 rng;
  pcg32_seed(&rng, 19, 29);

  Vec3 ref_point = vec3(4, 0.2, 0);
  float radius = 0.2f;

  for (int a = -11; a < 11; a++)
    for (int b = -11; b < 11; b++) {
      float choose_material = pcg32_f32(&rng);
      Vec3 center = vec3((float)a + 0.9f * pcg32_f32(&rng), radius, (float)b + 0.9f * pcg32_f32(&rng));

      if (vec3_length(vec3_sub(center, ref_point)) > 0.9f) {
        if (choose_material < 0.8f) {
          Vec3 color = vec3_mul(vec3_rand(&rng), vec3_rand(&rng));
          mat = Lambertian_new(Solid_new(color));
        } else if (choose_material < 0.95f) {
          Vec3 color = vec3_rand_between(&rng, 0.5f, 1);
          mat = Metal_new(Solid_new(color), pcg32_f32(&rng) * 0.5f);
        } else {
          mat = Dielectric_new(1.5f);
        }

        HittableList_append(&world->objects, Sphere_new(center, radius, mat));
      }
    }

  Hittable *bvh = BVHNode_new(&world->objects, &rng);
  free(world->objects.items);
  HittableList_init(&world->objects, 1);
  HittableList_append(&world->objects, bvh);

  camera->vfov = 20.0f;
  camera->background = vec3(0.7, 0.8, 1);
  camera->look_from = vec3(13, 2, 3);
  camera->look_to = VEC3_ZERO;
  camera->dof_angle = 0.6f;
}

void scene_checker(World *world, Camera *camera) {
  World_init(world, 2);

  Texture *checker_p = Checker_new(0.01f, Solid_new(vec3(0.2, 0.3, 0.1)), Solid_new(vec3(0.9, 0.9, 0.9)));
  Material *mat = Lambertian_new(checker_p);
  HittableList_append(&world->objects, Sphere_new(vec3(0, -10, 0), 10, mat));
  HittableList_append(&world->objects, Sphere_new(vec3(0, 10, 0), 10, mat));

  camera->vfov = 20.0f;
  camera->background = vec3(0.7, 0.8, 1);
  camera->look_from = vec3(13, 2, 3);
  camera->look_to = VEC3_ZERO;
}

void scene_earth(World *world, Camera *camera) {
  World_init(world, 1);

  Material *mat = Lambertian_new(Image_new("earthmap.jpg"));
  HittableList_append(&world->objects, Sphere_new(VEC3_ZERO, 2, mat));

  camera->vfov = 20.0f;
  camera->background = vec3(0.7, 0.8, 1);
  camera->look_from = vec3(13, 2, 3);
  camera->look_to = VEC3_ZERO;
}

void scene_perlin(World *world, Camera *camera) {
  World_init(world, 2);

  PCG32 rng;
  pcg32_seed(&rng, 19, 29);

  Material *mat = Lambertian_new(Perlin_new(4.0f, 7, &rng));
  HittableList_append(&world->objects, Sphere_new(vec3(0, -1000, 0), 1000, mat));
  HittableList_append(&world->objects, Sphere_new(vec3(0, 2, 0), 2, mat));

  camera->vfov = 20.0f;
  camera->background = vec3(0.7, 0.8, 1);
  camera->look_from = vec3(13, 2, 3);
  camera->look_to = VEC3_ZERO;
  camera->dof_angle = 0.0f;
  camera->focal_length = 10.0f;
}

void scene_simple_light(World *world, Camera *camera) {
  World_init(world, 4);

  PCG32 rng;
  pcg32_seed(&rng, 19, 29);

  Material *perlin = Lambertian_new(Perlin_new(4.0f, 7, &rng));
  Material *light = DiffuseLight_new(Solid_new(vec3(4, 4, 4)));
  HittableList_append(&world->objects, Sphere_new(vec3(0, -1000, 0), 1000, perlin));
  HittableList_append(&world->objects, Sphere_new(vec3(0, 2, 0), 2, perlin));

  Hittable *light_src;
  light_src = Quad_new(vec3(3, 1, -2), vec3(2, 0, 0), vec3(0, 2, 0), light);
  HittableList_append(&world->objects, light_src);
  HittableList_append(&world->lights, light_src);

  light_src = Sphere_new(vec3(0, 7, 0), 2, light);
  HittableList_append(&world->objects, light_src);
  HittableList_append(&world->lights, light_src);

  camera->vfov = 20.0f;
  camera->background = VEC3_ZERO;
  camera->look_from = vec3(26, 3, 6);
  camera->look_to = vec3(0, 2, 0);
}

void scene_cornell_box(World *world, Camera *camera) {
  World_init(world, 6 + 2);

  Material *red = Lambertian_new(Solid_new(vec3(0.65, 0.05, 0.05)));
  Material *white = Lambertian_new(Solid_new(vec3(0.73, 0.73, 0.73)));
  Material *green = Lambertian_new(Solid_new(vec3(0.12, 0.45, 0.15)));
  Material *light = DiffuseLight_new(Solid_new(vec3(15, 15, 15)));

  HittableList_append(&world->objects, Quad_new(vec3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), green));
  HittableList_append(&world->objects, Quad_new(vec3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), red));
  HittableList_append(&world->objects, Quad_new(vec3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), white));
  HittableList_append(&world->objects, Quad_new(vec3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), white));
  HittableList_append(&world->objects, Quad_new(vec3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), white));

  Hittable *light_src = Quad_new(vec3(343, 554, 332), vec3(-130, 0, 0), vec3(0, 0, -105), light);
  HittableList_append(&world->objects, light_src);
  HittableList_append(&world->lights, light_src);

  Hittable *box1 = Box_new(vec3(0, 0, 0), vec3(165, 330, 165), white);
  box1 = RotateY_new(box1, 15);
  box1 = Translate_new(box1, vec3(265, 0, 295));
  HittableList_append(&world->objects, box1);

  Hittable *box2 = Box_new(vec3(0, 0, 0), vec3(165, 165, 165), white);
  box2 = RotateY_new(box2, -18);
  box2 = Translate_new(box2, vec3(130, 0, 65));
  HittableList_append(&world->objects, box2);

  camera->aspect_ratio = 1.0f;
  camera->background = vec3(0, 0, 0);
  camera->vfov = 40.0f;
  camera->look_from = vec3(278, 278, -800);
  camera->look_to = vec3(278, 278, 0);
}

void scene_book2_final(World *world, Camera *camera, bool enable_bvh) {
  World_init(world, 11);

  PCG32 rng;
  pcg32_seed(&rng, 19, 29);

  int boxes_per_side = 20;
  Material *ground = Lambertian_new(Solid_new(vec3(0.48, 0.83, 0.53)));
  HittableList *boxes1_list = (HittableList *)HittableList_new(boxes_per_side * boxes_per_side);
  for (int i = 0; i < boxes_per_side; i++)
    for (int j = 0; j < boxes_per_side; j++) {
      float w = 100.0f;
      Vec3 p0 = vec3(-1000.0f + i * w, 0.0f, -1000.0f + j * w);
      Vec3 p1 = vec3(-1000.0f + (i + 1) * w, pcg32_f32_between(&rng, 1, 101), -1000.0f + (j + 1) * w);
      HittableList_append(boxes1_list, Box_new(p0, p1, ground));
    }

  Hittable *boxes1;
  if (enable_bvh) {
    boxes1 = BVHNode_new(boxes1_list, &rng);
    free(boxes1_list->items);
    free(boxes1_list);
  } else
    boxes1 = (Hittable *)boxes1_list;
  HittableList_append(&world->objects, boxes1);

  Material *light = DiffuseLight_new(Solid_new(vec3(7, 7, 7)));
  Hittable *light_src = Quad_new(vec3(123, 554, 147), vec3(300, 0, 0), vec3(0, 0, 265), light);
  HittableList_append(&world->objects, light_src);
  HittableList_append(&world->lights, light_src);

  Vec3 center1 = vec3(400, 400, 200);
  // Vec3 center2 = {430, 400, 200}; // TODO: moving sphere
  Material *sphere_mat = Lambertian_new(Solid_new(vec3(0.7, 0.3, 0.1)));
  HittableList_append(&world->objects, Sphere_new(center1, 50, sphere_mat));

  Material *glass = Dielectric_new(1.5);
  HittableList_append(&world->objects, Sphere_new(vec3(260, 150, 45), 50, glass));

  Material *metal = Metal_new(Solid_new(vec3(0.8, 0.8, 0.9)), 1.0);
  HittableList_append(&world->objects, Sphere_new(vec3(0, 150, 145), 50, metal));

  // subsurface material
  Hittable *boundary = Sphere_new(vec3(360, 150, 145), 70, glass);
  HittableList_append(&world->objects, boundary);
  HittableList_append(&world->objects, ConstantMedium_new(boundary, 0.2, Solid_new(vec3(0.2, 0.4, 0.9))));

  // mist
  boundary = Sphere_new(vec3(0, 0, 0), 5000, glass);
  HittableList_append(&world->objects, ConstantMedium_new(boundary, 0.0001, Solid_new(vec3(1, 1, 1))));

  Material *earth = Lambertian_new(Image_new("earthmap.jpg"));
  HittableList_append(&world->objects, Sphere_new(vec3(400, 200, 400), 100, earth));

  Material *perlin = Lambertian_new(Perlin_new(0.1, 7, &rng));
  HittableList_append(&world->objects, Sphere_new(vec3(220, 280, 300), 80, perlin));

  int ns = 1000;
  Material *white = Lambertian_new(Solid_new(vec3(0.73, 0.73, 0.73)));
  HittableList *boxes2_list = (HittableList *)HittableList_new(ns);
  for (int i = 0; i < ns; i++) {
    Vec3 center = vec3_rand_between(&rng, 0, 165);
    HittableList_append(boxes2_list, Sphere_new(center, 10, white));
  }
  Hittable *boxes2;
  if (enable_bvh) {
    boxes2 = BVHNode_new(boxes2_list, &rng);
    free(boxes2_list->items);
    free(boxes2_list);
  } else
    boxes2 = (Hittable *)boxes2_list;

  boxes2 = RotateY_new(boxes2, 15.0f);
  boxes2 = Translate_new(boxes2, vec3(-100, 270, 395));
  HittableList_append(&world->objects, boxes2);

  camera->aspect_ratio = 1.0f;
  camera->background = VEC3_ZERO;
  camera->vfov = 40.0f;
  camera->look_from = vec3(478, 278, -600);
  camera->look_to = vec3(278, 278, 0);
}

int main(int argc, char *argv[]) {
  assert(argc > 1);

  World world = {0};
  Camera camera;
  camera.aspect_ratio = 16.0f / 9.0f;
  camera.img_width = 500;
  camera.samples_per_pixel = 100;
  camera.max_depth = 50;
  camera.vup = vec3(0, 1, 0);
  camera.dof_angle = 0.0f;
  camera.focal_length = 10.0f;
  camera.lights_sampling_prob = 0.5f;

  if (argc > 3)
    camera.img_width = strtol(argv[2], NULL, 10);
  if (argc > 4)
    camera.samples_per_pixel = strtol(argv[3], NULL, 10);

  switch (strtol(argv[1], NULL, 10)) {
  default:
    fprintf(stderr, "Unsupported option. Default to 0\n");
  case 0:
    fprintf(stderr, "Book 1: Metal and Lambertian\n");
    scene_metal_and_lambertian(&world, &camera);
    break;
  case 1:
    fprintf(stderr, "Book 1: Final scene\n");
    scene_book1_final(&world, &camera);
    break;
  case 2:
    fprintf(stderr, "Book 2: Checker\n");
    scene_checker(&world, &camera);
    break;
  case 3:
    fprintf(stderr, "Book 2: Earth\n");
    scene_earth(&world, &camera);
    break;
  case 4:
    fprintf(stderr, "Book 2: Perlin noise\n");
    scene_perlin(&world, &camera);
    break;
  case 5:
    fprintf(stderr, "Book 2: Simple light\n");
    scene_simple_light(&world, &camera);
    break;
  case 6:
    fprintf(stderr, "Book 2: Cornell box\n");
    scene_cornell_box(&world, &camera);
    break;
  case 7:
    fprintf(stderr, "Book 2: Final scene\n");
    scene_book2_final(&world, &camera, true);
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
