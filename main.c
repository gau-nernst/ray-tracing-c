#include "tiff.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int width = 256;
  int height = 256;

  Image8 image = {width, height, 3, malloc(width * height * 3)};
  if (image.data == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    return 1;
  }

  for (int j = 0; j < height; j++)
    for (int i = 0; i < width; i++) {
      float r = (float)i / (width - 1);
      float g = (float)j / (height - 1);
      float b = 0;

      image.data[(j * width + i) * 3] = (int)(255.999 * r);
      image.data[(j * width + i) * 3 + 1] = (int)(255.999 * g);
      image.data[(j * width + i) * 3 + 2] = (int)(255.999 * b);
    }

  FILE *f = fopen("output.tiff", "wb");
  if (f == NULL) {
    fprintf(stderr, "Failed to open file\n");
    return 1;
  }
  write_tiff(f, image);

  return 0;
}
