#include "texture.h"
#include "utils.h"
#include "vec3.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static Vec3 Solid_value(const Texture *self_, float u, float v, Vec3 p) { return ((Solid *)self_)->color; }
void Solid_init(Solid *self, Vec3 color) { *self = (Solid){{Solid_value}, color}; }
Texture *Solid_new(Vec3 color) define_init_new(Solid, color);

static Vec3 Checker_value(const Texture *self_, float u, float v, Vec3 p) {
  const Checker *self = (const Checker *)self_;
  // int x = (int)floorf(p.x / texture->scale);
  // int y = (int)floorf(p.y / texture->scale);
  // int z = (int)floorf(p.z / texture->scale);
  // return texture_value((x + y + z) % 2 ? texture->odd : texture->even, u, v, p);
  int int_u = (int)floorf(u / self->scale);
  int int_v = (int)floorf(v / self->scale);
  Texture *select = (int_u + int_v) % 2 ? self->odd : self->even;
  return select->value(select, u, v, p);
}
void Checker_init(Checker *self, float scale, Texture *even, Texture *odd) {
  *self = (Checker){{Checker_value}, scale, even, odd};
}
Texture *Checker_new(float scale, Texture *even, Texture *odd) define_init_new(Checker, scale, even, odd);

static Vec3 Image_value(const Texture *self_, float u, float v, Vec3 p) {
  const Image *self = (const Image *)self_;
  // nearest neighbour sampling
  // flip v to image coordinates
  int i = (int)roundf(u * (float)(self->width - 1));
  int j = (int)roundf((1.0f - v) * (float)(self->height - 1));
  int offset = ((j * self->width) + i) * 3;
  return vec3((float)self->buffer[offset] / 255.0f, (float)self->buffer[offset + 1] / 255.0f,
              (float)self->buffer[offset + 2] / 255.0f);
}
void Image_init(Image *image, char *filename) {
  image->texture.value = Image_value;
  image->buffer = stbi_load(filename, &image->width, &image->height, NULL, 3);
  assert((image->buffer != NULL) && "Unable to read image");
}
Texture *Image_new(char *filename) define_init_new(Image, filename);

static void Perlin_permute(int perm[N_PERLIN], PCG32 *rng);
static float Perlin_turbulence(const Perlin *perlin, Vec3 p);
static Vec3 Perlin_value(const Texture *self_, float u, float v, Vec3 p) {
  const Perlin *self = (const Perlin *)self_;
  p = vec3_mul(p, self->scale);
  float value = 0.5f * (1.0f + sinf(p.z + 10.0f * Perlin_turbulence(self, p)));
  return vec3(value, value, value);
}
void Perlin_init(Perlin *perlin, float scale, int depth, PCG32 *rng) {
  perlin->texture.value = Perlin_value;
  perlin->scale = scale;
  perlin->depth = depth;
  for (int i = 0; i < N_PERLIN; i++)
    perlin->grad_field[i] = vec3_rand_unit_vector(rng);
  Perlin_permute(perlin->perm_x, rng);
  Perlin_permute(perlin->perm_y, rng);
  Perlin_permute(perlin->perm_z, rng);
}
Texture *Perlin_new(float scale, int depth, PCG32 *rng) define_init_new(Perlin, scale, depth, rng);

static void Perlin_permute(int perm[N_PERLIN], PCG32 *rng) {
  for (int i = 0; i < N_PERLIN; i++)
    perm[i] = i;
  for (int i = N_PERLIN - 1; i > 0; i--) {
    uint32_t target = pcg32_u32_between(rng, 0, i + 1);
    int tmp = perm[i];
    perm[i] = perm[target];
    perm[target] = tmp;
  }
}

static float hermitian_smoothing(float t) { return t * t * (3.0f - 2.0f * t); }

static float Perlin_noise(const Perlin *perlin, Vec3 p) {
  int i = (int)floorf(p.x);
  int j = (int)floorf(p.y);
  int k = (int)floorf(p.z);

  float t1 = p.x - (float)i;
  float t2 = p.y - (float)j;
  float t3 = p.z - (float)k;

  float tt1 = hermitian_smoothing(t1);
  float tt2 = hermitian_smoothing(t2);
  float tt3 = hermitian_smoothing(t3);

  float value = 0;
  for (int di = 0; di < 2; di++)
    for (int dj = 0; dj < 2; dj++)
      for (int dk = 0; dk < 2; dk++) {
        Vec3 grad = perlin->grad_field[perlin->perm_x[(i + di) & 255] ^ perlin->perm_y[(j + dj) & 255] ^
                                       perlin->perm_z[(k + dk) & 255]];
        Vec3 weight = vec3(t1 - di, t2 - dj, t3 - dk);
        value += vec3_dot(grad, weight) * (di * tt1 + (1 - di) * (1.0f - tt1)) * (dj * tt2 + (1 - dj) * (1.0f - tt2)) *
                 (dk * tt3 + (1 - dk) * (1.0f - tt3));
      }

  return value;
}

static float Perlin_turbulence(const Perlin *perlin, Vec3 p) {
  float value = 0.0f;
  float weight = 1.0f;
  for (int i = 0; i < perlin->depth; i++) {
    value += weight * Perlin_noise(perlin, p);
    weight *= 0.5f;
    p = vec3_mul(p, 2.0f);
  }
  return fabsf(value);
}
