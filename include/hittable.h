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

typedef union AABB {
  struct {
    float x[2];
    float y[2];
    float z[2];
  };
  float values[3][2];
} AABB;

typedef struct Hittable Hittable;
typedef struct HittableVTable {
  bool (*hit)(const Hittable *self, const Ray *ray, float t_min, float t_max, HitRecord *rec, PCG32 *rng);
} HittableVTable;

struct Hittable {
  HittableVTable *vtable;
  AABB bbox;
};

typedef struct HittableList {
  Hittable hittable;
  size_t max_size;
  size_t size;
  Hittable **items;
} HittableList;

void HittableList_init(HittableList *self, size_t max_size);
Hittable *HittableList_new(size_t max_size);
void HittableList_append(HittableList *self, Hittable *item);

typedef struct Sphere {
  Hittable hittable;
  Vec3 center;
  float radius;
  Material *material;
} Sphere;

void Sphere_init(Sphere *self, Vec3 center, float radius, Material *mat);
Hittable *Sphere_new(Vec3 center, float radius, Material *mat);

typedef struct Quad {
  Hittable hittable;
  Vec3 Q;
  Vec3 u;
  Vec3 v;
  Vec3 normal;
  float D;
  Vec3 w;
  Material *material;
} Quad;

void Quad_init(Quad *self, Vec3 Q, Vec3 u, Vec3 v, Material *mat);
Hittable *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material *mat);
Hittable *Box_new(Vec3 a, Vec3 b, Material *mat);

typedef struct BVHNode {
  Hittable hittable;
  Hittable *left;
  Hittable *right;
} BVHNode;

void BVHNode_init(BVHNode *self, const HittableList *list, PCG32 *rng);
Hittable *BVHNode_new(const HittableList *list, PCG32 *rng);

typedef struct Translate {
  Hittable hittable;
  Hittable *object;
  Vec3 offset;
} Translate;

void Translate_init(Translate *self, Hittable *object, Vec3 offset);
Hittable *Translate_new(Hittable *object, Vec3 offset);

typedef struct RotateY {
  Hittable hittable;
  Hittable *object;
  float sin_theta;
  float cos_theta;
} RotateY;

void RotateY_init(RotateY *self, Hittable *object, float angle);
Hittable *RotateY_new(Hittable *object, float angle);

typedef struct ConstantMedium {
  Hittable hittable;
  Hittable *boundary;
  float neg_inv_density;
  Material *phase_fn;
} ConstantMedium;

void ConstantMedium_init(ConstantMedium *self, Hittable *boundary, float density, Texture *albedo);
Hittable *ConstantMedium_new(Hittable *boundary, float density, Texture *albedo);

#endif // HITTABLE_H
