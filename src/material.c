#include "material.h"
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

void texture_image_load(Texture *texture, char *filename) {
  texture->type = IMAGE;
  texture->image = stbi_load(filename, &texture->width, &texture->height, NULL, 3);
}

Vec3 texture_value_checker(Texture *texture, float u, float v, Vec3 p) {
  // int x = (int)floorf(p.x / texture->scale);
  // int y = (int)floorf(p.y / texture->scale);
  // int z = (int)floorf(p.z / texture->scale);
  // return texture_value((x + y + z) % 2 ? texture->odd : texture->even, u, v, p);
  int int_u = (int)floorf(u / texture->scale);
  int int_v = (int)floorf(v / texture->scale);
  return texture_value((int_u + int_v) % 2 ? texture->odd : texture->even, u, v, p);
}

Vec3 texture_value_image(Texture *texture, float u, float v, Vec3 p) {
  // nearest neighbour sampling
  // flip v to image coordinates
  int i = (int)roundf(u * (float)(texture->width - 1));
  int j = (int)roundf((1.0f - v) * (float)(texture->height - 1));
  return (Vec3){
      (float)texture->image[((j * texture->width) + i) * 3] / 255.0f,
      (float)texture->image[((j * texture->width) + i) * 3 + 1] / 255.0f,
      (float)texture->image[((j * texture->width) + i) * 3 + 2] / 255.0f,
  };
}

Vec3 texture_value(Texture *texture, float u, float v, Vec3 p) {
  switch (texture->type) {
  case SOLID:
    return texture->color;
  case CHECKER:
    return texture_value_checker(texture, u, v, p);
  case IMAGE:
    return texture_value_image(texture, u, v, p);
  default:
    return (Vec3){0.0f, 0.0f, 0.0f};
  }
}

bool scatter_normal(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = vec3float_mul(vec3_add(hit_record->normal, 1.0f), 0.5f); // [-1,1] -> [0,1]
  return false;
}

bool scatter_lambertian(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
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
  if (eta * sin_theta > 1.0f || schlick_r > pcg32_randomf_r(rng)) {
    *scattered = reflect(incident, hit_record->normal);
  } else {
    Vec3 r_perp = vec3_mul(vec3_add(incident, vec3_mul(hit_record->normal, cos_theta)), eta);
    Vec3 r_para = vec3_mul(hit_record->normal, -sqrtf(fabsf(1.0f - vec3_length2(r_perp))));
    *scattered = vec3_add(r_perp, r_para);
  }
  *color = texture_value(hit_record->material->albedo, hit_record->u, hit_record->v, hit_record->p);
  return true;
}

bool scatter_default(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  *color = (Vec3){0.0f, 0.0f, 0.0f};
  return false;
}

bool scatter(Vec3 incident, HitRecord *hit_record, PCG32State *rng, Vec3 *scattered, Vec3 *color) {
  switch (hit_record->material->type) {
  case NORMAL:
    return scatter_normal(incident, hit_record, rng, scattered, color);
  case LAMBERTIAN:
    return scatter_lambertian(incident, hit_record, rng, scattered, color);
  case METAL:
    return scatter_metal(incident, hit_record, rng, scattered, color);
  case DIELECTRIC:
    return scatter_dielectric(incident, hit_record, rng, scattered, color);
  default:
    return scatter_default(incident, hit_record, rng, scattered, color);
  }
}
