#include "pcg32.h"
#include "vec3.h"
#include <stdint.h>

typedef struct Texture {
  enum TextureType {
    SOLID,
    CHECKER,
    IMAGE,
    PERLIN,
  } type;
  void *ptr;
} Texture;

Vec3 texture_value(Texture texture, float u, float v, Vec3 p);

typedef struct Checker {
  float scale;
  Texture even;
  Texture odd;
} Checker;

typedef struct Image {
  int width;
  int height;
  uint8_t *buffer;
} Image;

void image_load(Image *image, char *filename);

#define N_PERLIN 256

typedef struct Perlin {
  float scale;
  int depth;
  Vec3 grad_field[N_PERLIN];
  int perm_x[N_PERLIN];
  int perm_y[N_PERLIN];
  int perm_z[N_PERLIN];
} Perlin;

void perlin_init(Perlin *perlin, PCG32State *rng);
