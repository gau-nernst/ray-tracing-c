#ifndef TIFF_H
#define TIFF_H

#include <stdint.h>
#include <stdio.h>

int write_tiff(FILE *f, int width, int height, int n_channels, uint8_t *buffer);

#endif // TIFF_H
