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

float vec3_length2(Vec3 u) { return vec3_dot(u, u); }
float vec3_length(Vec3 u) { return sqrtf(vec3_length2(u)); }
float vec3_dot(Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
Vec3 vec3_cross(Vec3 u, Vec3 v) { return (Vec3){u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y}; }
Vec3 vec3_unit(Vec3 u) { return vec3float_mul(u, 1.0f / vec3_length(u)); }

Vec3 ray_at(Ray ray, float t) { return vec3vec3_add(ray.origin, vec3float_mul(ray.direction, t)); }

HitRecord hit_sphere(const Sphere *sphere, const Ray *ray, float t_min, float t_max) {
  HitRecord hit_record = {0};

  Vec3 oc = vec3_sub(ray->origin, sphere->center);
  float a = vec3_length2(ray->direction);
  float b = vec3_dot(oc, ray->direction);
  float c = vec3_length2(oc) - sphere->radius * sphere->radius;
  float disc = b * b - a * c;

  if (disc < 0)
    return hit_record;

  float disc_sqrt = sqrtf(disc);
  float root = (-b - disc_sqrt) / a;
  if (root <= t_min || root >= t_max) {
    root = (-b + disc_sqrt) / a;
    if (root <= t_min || root >= t_max)
      return hit_record;
  }

  hit_record.t = root;
  hit_record.p = ray_at(*ray, root);
  hit_record.normal = vec3_mul(vec3_sub(hit_record.p, sphere->center), 1.0f / sphere->radius);
  return hit_record;
}
