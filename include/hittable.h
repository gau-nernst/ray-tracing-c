#ifndef HITTABLE_H
#define HITTABLE_H

#include "material.h"
#include "vec3.h"
#include <stdbool.h>

typedef struct Ray {
  Vec3 origin;
  Vec3 direction;
} Ray;

Vec3 ray_at(const Ray *ray, float t);

typedef struct AABB {
  float x[3][2];
} AABB;

AABB AABB_from_Vec3(Vec3 a, Vec3 b);
AABB AABB_from_AABB(AABB *a, AABB *b);

typedef enum HittableType {
  HITTABLE_LIST,
  SPHERE,
  QUAD,
  BVH_NODE,
  TRANSLATE,
  ROTATE_Y,
  CONSTANT_MEDIUM,
} HittableType;

typedef struct Hittable {
  HittableType type;
  void *ptr;
} Hittable;

#define hittable(ptr)                                                                                                  \
  (Hittable) {                                                                                                         \
    _Generic((ptr),                                                                                                    \
        HittableList *: HITTABLE_LIST,                                                                                 \
        Sphere *: SPHERE,                                                                                              \
        Quad *: QUAD,                                                                                                  \
        BVHNode *: BVH_NODE,                                                                                           \
        Translate *: TRANSLATE,                                                                                        \
        RotateY *: ROTATE_Y,                                                                                           \
        ConstantMedium *: CONSTANT_MEDIUM),                                                                            \
        ptr                                                                                                            \
  }

typedef struct HittableList {
  size_t max_size;
  size_t size;
  Hittable *items;
  AABB bbox;
} HittableList;

void HittableList_init(HittableList *list, size_t max_size);
HittableList *HittableList_new(size_t max_size);
void HittableList_append(HittableList *list, Hittable item);

bool Hittable_hit(Hittable obj, const Ray *ray, float t_min, float t_max, HitRecord *hit_record, PCG32State *rng);
bool HittableList_hit(const HittableList *list, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                      PCG32State *rng);

typedef struct Sphere {
  Vec3 center;
  float radius;
  Material material;
  AABB bbox;
} Sphere;

void Sphere_init(Sphere *sphere, Vec3 center, float radius, Material material);
Sphere *Sphere_new(Vec3 center, float radius, Material material);

typedef struct Quad {
  Vec3 Q;
  Vec3 u;
  Vec3 v;
  Vec3 normal;
  float D;
  Vec3 w;
  Material material;
} Quad;

void Quad_init(Quad *quad, Vec3 Q, Vec3 u, Vec3 v, Material mat);
Quad *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material mat);

HittableList *Box_new(Vec3 a, Vec3 b, Material mat);

typedef struct BVHNode {
  Hittable left;
  Hittable right;
  AABB bbox;
} BVHNode;

typedef struct Translate {
  Hittable object;
  Vec3 offset;
} Translate;

Translate *Translate_new(Hittable object, Vec3 offset);

typedef struct RotateY {
  Hittable object;
  float sin_theta;
  float cos_theta;
} RotateY;

void RotateY_init(RotateY *rotate_y, Hittable object, float angle);
RotateY *RotateY_new(Hittable object, float angle);

typedef struct ConstantMedium {
  Hittable boundary;
  float neg_inv_density;
  Material phase_fn;
} ConstantMedium;

void ConstantMedium_init(ConstantMedium *constant_medium, Hittable boundary, float density, Texture albedo);
ConstantMedium *ConstantMedium_new(Hittable boundary, float density, Texture albedo);

#endif // HITTABLE_H
