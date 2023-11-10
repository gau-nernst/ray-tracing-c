#include "material.h"
#include <math.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

Vec3 checker_texture_value(Checker *checker, float u, float v, Vec3 p) {
  // int x = (int)floorf(p.x / texture->scale);
  // int y = (int)floorf(p.y / texture->scale);
  // int z = (int)floorf(p.z / texture->scale);
  // return texture_value((x + y + z) % 2 ? texture->odd : texture->even, u, v, p);
  int int_u = (int)floorf(u / checker->scale);
  int int_v = (int)floorf(v / checker->scale);
  return texture_value((int_u + int_v) % 2 ? checker->odd : checker->even, u, v, p);
}

void image_load(Image *image, char *filename) {
  image->buffer = stbi_load(filename, &image->width, &image->height, NULL, 3);
}

Vec3 image_texture_value(Image *image, float u, float v) {
  // nearest neighbour sampling
  // flip v to image coordinates
  int i = (int)roundf(u * (float)(image->width - 1));
  int j = (int)roundf((1.0f - v) * (float)(image->height - 1));
  int offset = ((j * image->width) + i) * 3;
  return (Vec3){
      (float)image->buffer[offset] / 255.0f,
      (float)image->buffer[offset + 1] / 255.0f,
      (float)image->buffer[offset + 2] / 255.0f,
  };
}

void perlin_permute(int perm[N_PERLIN], PCG32State *rng) {
  for (int i = 0; i < N_PERLIN; i++)
    perm[i] = i;
  for (int i = N_PERLIN - 1; i > 0; i--) {
    uint32_t target = pcg32_u32_between(rng, 0, i + 1);
    int tmp = perm[i];
    perm[i] = perm[target];
    perm[target] = tmp;
  }
}

void perlin_init(Perlin *perlin, PCG32State *rng) {
  for (int i = 0; i < N_PERLIN; i++)
    perlin->grad_field[i] = vec3_rand_unit_vector(rng);
  perlin_permute(perlin->perm_x, rng);
  perlin_permute(perlin->perm_y, rng);
  perlin_permute(perlin->perm_z, rng);
}

float hermitian_smoothing(float t) { return t * t * (3.0f - 2.0f * t); }

float perlin_noise(Perlin *perlin, Vec3 p) {
  int i = (int)floorf(p.x[0]);
  int j = (int)floorf(p.x[1]);
  int k = (int)floorf(p.x[2]);

  float t1 = p.x[0] - (float)i;
  float t2 = p.x[1] - (float)j;
  float t3 = p.x[2] - (float)k;

  float tt1 = hermitian_smoothing(t1);
  float tt2 = hermitian_smoothing(t2);
  float tt3 = hermitian_smoothing(t3);

  float value = 0;
  for (int di = 0; di < 2; di++)
    for (int dj = 0; dj < 2; dj++)
      for (int dk = 0; dk < 2; dk++) {
        Vec3 grad = perlin->grad_field[perlin->perm_x[(i + di) & 255] ^ perlin->perm_y[(j + dj) & 255] ^
                                       perlin->perm_z[(k + dk) & 255]];
        Vec3 weight = {t1 - di, t2 - dj, t3 - dk};
        value += vec3_dot(grad, weight) * (di * tt1 + (1 - di) * (1.0f - tt1)) * (dj * tt2 + (1 - dj) * (1.0f - tt2)) *
                 (dk * tt3 + (1 - dk) * (1.0f - tt3));
      }

  return value;
}

float perlin_turbulence(Perlin *perlin, Vec3 p) {
  float value = 0.0f;
  float weight = 1.0f;
  for (int i = 0; i < perlin->depth; i++) {
    value += weight * perlin_noise(perlin, p);
    weight *= 0.5f;
    p = vec3_mul(p, 2.0f);
  }
  return fabs(value);
}

Vec3 perlin_texture_value(Perlin *perlin, Vec3 p) {
  p = vec3_mul(p, perlin->scale);
  // float value = perlin_noise(perlin, p);
  // value = (value + 1.0f) * 0.5f;
  float value = perlin_turbulence(perlin, p);
  value = (sinf(10.0f * value + p.x[2]) + 1.0f) * 0.5f;
  return (Vec3){value, value, value};
}

Vec3 texture_value(Texture texture, float u, float v, Vec3 p) {
  switch (texture.type) {
  case SOLID:
    return *texture.color;
  case CHECKER:
    return checker_texture_value(texture.checker, u, v, p);
  case IMAGE:
    return image_texture_value(texture.image, u, v);
  case PERLIN:
    return perlin_texture_value(texture.perlin, p);
  default:
    return (Vec3){0.0f, 0.0f, 0.0f};
  }
}

bool scatter_lambertian(HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  Vec3 new_direction = vec3_add(hit_record->normal, vec3_rand_unit_vector(rng));
  if (vec3_near_zero(new_direction))
    new_direction = hit_record->normal;
  *scattered = new_direction;
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
  return vec3_sub(incident, vec3_mul(normal, 2.0f * vec3_dot(incident, normal)));
}

bool scatter_metal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *scattered =
      vec3_add(reflect(incident, hit_record->normal), vec3_mul(vec3_rand_unit_vector(rng), hit_record->material->fuzz));
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return vec3_dot(*scattered, hit_record->normal) > 0; // check for degeneration
}

bool scatter_dielectric(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  float eta = hit_record->front_face ? 1.0f / hit_record->material->eta : hit_record->material->eta;
  incident = vec3_unit(incident);

  float cos_theta = fminf(-vec3_dot(incident, hit_record->normal), 1.0f);
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

  float schlick_r = (1.0f - eta) / (1.0f + eta);
  schlick_r *= schlick_r;
  schlick_r += (1 - schlick_r) * powf(1.0f - cos_theta, 5.0f);

  // total internal reflection or partial reflections (Fresnel)
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_f32(rng)) {
    *scattered = reflect(incident, hit_record->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(incident, vec3_mul(hit_record->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(hit_record->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3_add(r_perp, r_para);
  }
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  switch (hit_record->material->type) {
  case NORMAL:
    *color = vec3float_mul(vec3_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
    return false;
  case LAMBERTIAN:
    return scatter_lambertian(hit_record, rng, scattered, color);
  case METAL:
    return scatter_metal(incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return scatter_dielectric(incident, hit_record, rng, scattered, color);
  default:
    *color = (Vec3){0.0f, 0.0f, 0.0f};
    return false;
  }
}

Vec3 emit(HitRecord *hit_record) {
  switch (hit_record->material->type) {
  case DIFFUSE_LIGHT:
    return *hit_record->material->albedo.color;
  default:
    return (Vec3){0.0f, 0.0f, 0.0f};
  }
}
