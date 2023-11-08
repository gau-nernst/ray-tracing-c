#ifndef AABB_H
#define AABB_H

#include "raytracing.h"
#include <stdbool.h>

typedef struct AABB {
  float x[3][2];
} AABB;

bool aabb_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max);

#endif // AABB_H
