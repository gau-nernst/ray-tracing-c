#include "pcg32.h"
#include "vec3.h"
#include <stdint.h>

typedef enum TextureType {
  SOLID,
  CHECKER,
  IMAGE,
  PERLIN,
} TextureType;

typedef struct Texture {
  TextureType type;
  void *ptr;
} Texture;

#define texture(ptr)                                                                                                   \
  (Texture) { _Generic((ptr), Vec3 *: SOLID, Checker *: CHECKER, Image *: IMAGE, Perlin *: PERLIN), ptr }

Vec3 texture_value(Texture texture, float u, float v, Vec3 p);

typedef struct Checker {
  float scale;
  Texture even;
  Texture odd;
} Checker;

Checker *Checker_new(float scale, Texture even, Texture odd);

typedef struct Image {
  int width;
  int height;
  uint8_t *buffer;
} Image;

void Image_init(Image *image, char *filename);
Image *Image_new(char *filename);

#define N_PERLIN 256

typedef struct Perlin {
  float scale;
  int depth;
  Vec3 grad_field[N_PERLIN];
  int perm_x[N_PERLIN];
  int perm_y[N_PERLIN];
  int perm_z[N_PERLIN];
} Perlin;

void Perlin_init(Perlin *perlin, float scale, int depth, PCG32State *rng);
Perlin *Perlin_new(float scale, int depth, PCG32State *rng);
