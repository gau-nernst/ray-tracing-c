#include "aabb.h"

bool aabb_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max) {
  for (int a = 0; a < 3; a++) {
    float invD = 1.0f / ray->direction.x[a];
    float t0 = (min(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;
    float t1 = (max(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;

    if (invD < 0) {
      float tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    if (t0 > t_min)
      t_min = t0;
    if (t1 < t_max)
      t_max = t1;
    if (t_max < t_min)
      return false;
  }
  return true;
}

AABB aabb_pad(const AABB *aabb) {
  float delta = 1e-4f;
  AABB padded;
  for (int a = 0; a < 3; a++) {
    if (aabb->x[a][1] - aabb->x[a][0] < delta) {
      padded.x[a][0] = aabb->x[a][0] - delta;
      padded.x[a][1] = aabb->x[a][1] + delta;
    } else {
      padded.x[a][0] = aabb->x[a][0];
      padded.x[a][1] = aabb->x[a][1];
    }
  }
  return padded;
}
