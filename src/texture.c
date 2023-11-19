#include "texture.h"
#include "utils.h"

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

Checker *Checker_new(float scale, Texture even, Texture odd) define_struct_new(Checker, scale, even, odd);
Image *Image_new(char *filename) define_init_new(Image, filename);
Perlin *Perlin_new(float scale, int depth, PCG32State *rng) define_init_new(Perlin, scale, depth, rng);

static Vec3 Solid_value(Vec3 *color);
static Vec3 Checker_value(Checker *checker, float u, float v, Vec3 p);
static Vec3 Image_value(Image *image, float u, float v);
static Vec3 Perlin_value(Perlin *perlin, Vec3 p);

Vec3 Texture_value(Texture texture, float u, float v, Vec3 p) {
  switch (texture.type) {
  case SOLID:
    return Solid_value(texture.ptr);
  case CHECKER:
    return Checker_value(texture.ptr, u, v, p);
  case IMAGE:
    return Image_value(texture.ptr, u, v);
  case PERLIN:
    return Perlin_value(texture.ptr, p);
  default:
    return (Vec3){0, 0, 0};
  }
}

static Vec3 Solid_value(Vec3 *color) { return *color; }

static Vec3 Checker_value(Checker *checker, float u, float v, Vec3 p) {
  // int x = (int)floorf(p.x / texture->scale);
  // int y = (int)floorf(p.y / texture->scale);
  // int z = (int)floorf(p.z / texture->scale);
  // return texture_value((x + y + z) % 2 ? texture->odd : texture->even, u, v, p);
  int int_u = (int)floorf(u / checker->scale);
  int int_v = (int)floorf(v / checker->scale);
  return Texture_value((int_u + int_v) % 2 ? checker->odd : checker->even, u, v, p);
}

void Image_init(Image *image, char *filename) {
  image->buffer = stbi_load(filename, &image->width, &image->height, NULL, 3);
}

static Vec3 Image_value(Image *image, float u, float v) {
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

static void Perlin_permute(int perm[N_PERLIN], PCG32State *rng) {
  for (int i = 0; i < N_PERLIN; i++)
    perm[i] = i;
  for (int i = N_PERLIN - 1; i > 0; i--) {
    uint32_t target = pcg32_u32_between(rng, 0, i + 1);
    int tmp = perm[i];
    perm[i] = perm[target];
    perm[target] = tmp;
  }
}

void Perlin_init(Perlin *perlin, float scale, int depth, PCG32State *rng) {
  perlin->scale = scale;
  perlin->depth = depth;
  for (int i = 0; i < N_PERLIN; i++)
    perlin->grad_field[i] = vec3_rand_unit_vector(rng);
  Perlin_permute(perlin->perm_x, rng);
  Perlin_permute(perlin->perm_y, rng);
  Perlin_permute(perlin->perm_z, rng);
}

static float hermitian_smoothing(float t) { return t * t * (3.0f - 2.0f * t); }

static float Perlin_noise(Perlin *perlin, Vec3 p) {
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

static float Perlin_turbulence(Perlin *perlin, Vec3 p) {
  float value = 0.0f;
  float weight = 1.0f;
  for (int i = 0; i < perlin->depth; i++) {
    value += weight * Perlin_noise(perlin, p);
    weight *= 0.5f;
    p = vec3_mul(p, 2.0f);
  }
  return fabs(value);
}

static Vec3 Perlin_value(Perlin *perlin, Vec3 p) {
  p = vec3_mul(p, perlin->scale);
  float value = 0.5f * (1.0f + sinf(p.x[2] + 10.0f * Perlin_turbulence(perlin, p)));
  return (Vec3){value, value, value};
}
