#include <stdint.h>
#include <stdio.h>

typedef struct Image8 {
  int width;
  int height;
  int n_channels;
  uint8_t *data;
} Image8;

int write_tiff(FILE *f, Image8 image);
