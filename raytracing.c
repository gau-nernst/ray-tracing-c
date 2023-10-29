#include "raytracing.h"
#include <math.h>

Vec3 *vec3vec3_add_(Vec3 *self, Vec3 other) {
  self->x += other.x;
  self->y += other.y;
  self->z += other.z;
  return self;
}

Vec3 *vec3vec3_sub_(Vec3 *self, Vec3 other) {
  self->x -= other.x;
  self->y -= other.y;
  self->z -= other.z;
  return self;
}

Vec3 *vec3float_mul_(Vec3 *self, float other) {
  self->x *= other;
  self->y *= other;
  self->z *= other;
  return self;
}

Vec3 vec3vec3_add(Vec3 u, Vec3 v) { return (Vec3){u.x + v.x, u.y + v.y, u.z + v.z}; }
Vec3 vec3vec3_sub(Vec3 u, Vec3 v) { return (Vec3){u.x - v.x, u.y - v.y, u.z - v.z}; }
Vec3 vec3vec3_mul(Vec3 u, Vec3 v) { return (Vec3){u.x * v.x, u.y * v.y, u.z * v.z}; }
Vec3 vec3float_add(Vec3 u, float v) { return (Vec3){u.x + v, u.y + v, u.z + v}; }
Vec3 vec3float_sub(Vec3 u, float v) { return (Vec3){u.x - v, u.y - v, u.z - v}; }
Vec3 vec3float_mul(Vec3 u, float v) { return (Vec3){u.x * v, u.y * v, u.z * v}; }

float vec3_length2(Vec3 u) { return u.x * u.x + u.y * u.y + u.z * u.z; }
float vec3_length(Vec3 u) { return sqrtf(vec3_length2(u)); }
float vec3_dot(Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y + u.z + v.z; }
Vec3 vec3_cross(Vec3 u, Vec3 v) { return (Vec3){u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y}; }
Vec3 vec3_unit(Vec3 u) { return vec3float_mul(u, 1.0f / vec3_length(u)); }

Vec3 ray_at(Ray ray, float t) { return vec3vec3_add(ray.origin, vec3float_mul(ray.direction, t)); }
