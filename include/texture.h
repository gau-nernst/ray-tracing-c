#include "pcg32.h"
#include "vec3.h"
#include <stdint.h>

typedef struct Texture Texture;
struct Texture {
  Vec3 (*value)(Texture *self, float u, float v, Vec3 p);
};

typedef struct Solid {
  Texture texture;
  Vec3 color;
} Solid;

void Solid_init(Solid *self, Vec3 color);
Texture *Solid_new(Vec3 color);

typedef struct Checker {
  Texture texture;
  float scale;
  Texture *even;
  Texture *odd;
} Checker;

void Checker_init(Checker *self, float scale, Texture *even, Texture *odd);
Texture *Checker_new(float scale, Texture *even, Texture *odd);

typedef struct Image {
  Texture texture;
  int width;
  int height;
  uint8_t *buffer;
} Image;

void Image_init(Image *self, char *filename);
Texture *Image_new(char *filename);

#define N_PERLIN 256

typedef struct Perlin {
  Texture texture;
  float scale;
  int depth;
  Vec3 grad_field[N_PERLIN];
  int perm_x[N_PERLIN];
  int perm_y[N_PERLIN];
  int perm_z[N_PERLIN];
} Perlin;

void Perlin_init(Perlin *self, float scale, int depth, PCG32 *rng);
Texture *Perlin_new(float scale, int depth, PCG32 *rng);
