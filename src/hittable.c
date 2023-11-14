#include "hittable.h"

#define _USE_MATH_DEFINES // for MSVC
#include <math.h>

Vec3 ray_at(const Ray *ray, float t) { return vec3vec3_add(ray->origin, vec3float_mul(ray->direction, t)); }

define_list_source(Hittable);

Sphere *Sphere_new(Vec3 center, float radius, Material mat) define_struct_new(Sphere, center, radius, mat);

void Quad_init(Quad *quad, Vec3 Q, Vec3 u, Vec3 v, Material material) {
  quad->Q = Q;
  quad->u = u;
  quad->v = v;
  quad->material = material;

  Vec3 n = vec3_cross(quad->u, quad->v);
  quad->normal = vec3_unit(n);
  quad->D = vec3_dot(quad->normal, quad->Q);
  quad->w = vec3_div(n, vec3_length2(n));
}
Quad *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material mat) define_init_new(Quad, Q, u, v, mat);

HittableList *Box_new(Vec3 a, Vec3 b, Material mat) {
  HittableList *box = HittableList_new(6);

  Vec3 min_p = vec3_min(a, b);
  Vec3 max_p = vec3_max(a, b);

  Vec3 dx = {max_p.x[0] - min_p.x[0], 0, 0};
  Vec3 dy = {0, max_p.x[1] - min_p.x[1], 0};
  Vec3 dz = {0, 0, max_p.x[2] - min_p.x[2]};

  HittableList_append(box, hittable(Quad_new((Vec3){min_p.x[0], min_p.x[1], max_p.x[2]}, dx, dy, mat))); // front
  HittableList_append(box,
                      hittable(Quad_new((Vec3){max_p.x[0], min_p.x[1], max_p.x[2]}, vec3_neg(dz), dy, mat))); // right
  HittableList_append(box,
                      hittable(Quad_new((Vec3){max_p.x[0], min_p.x[1], min_p.x[2]}, vec3_neg(dx), dy, mat))); // back
  HittableList_append(box, hittable(Quad_new((Vec3){min_p.x[0], min_p.x[1], min_p.x[2]}, dz, dy, mat)));      // left
  HittableList_append(box,
                      hittable(Quad_new((Vec3){min_p.x[0], max_p.x[1], max_p.x[2]}, dx, vec3_neg(dz), mat))); // top
  HittableList_append(box, hittable(Quad_new((Vec3){min_p.x[0], min_p.x[1], min_p.x[2]}, dx, dz, mat)));      // bottom

  return box;
}

Translate *Translate_new(Hittable object, Vec3 offset) define_struct_new(Translate, object, offset);

void RotateY_init(RotateY *rotate_y, Hittable object, float angle) {
  angle = angle * (float)M_PI / 180.f;
  *rotate_y = (RotateY){object, sinf(angle), cosf(angle)};
}
RotateY *RotateY_new(Hittable object, float angle) define_init_new(RotateY, object, angle);

void ConstantMedium_init(ConstantMedium *constant_medium, Hittable boundary, float density, Texture albedo) {
  *constant_medium = (ConstantMedium){
      boundary,
      -1.0f / density,
      material(Isotropic_new(albedo)),
  };
}
ConstantMedium *ConstantMedium_new(Hittable boundary, float density, Texture albedo)
    define_init_new(ConstantMedium, boundary, density, albedo);

// static bool aabb_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max) {
//   for (int a = 0; a < 3; a++) {
//     float invD = 1.0f / ray->direction.x[a];
//     float t0 = (min(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;
//     float t1 = (max(aabb->x[a][0], aabb->x[a][1]) - ray->origin.x[a]) * invD;

//     if (invD < 0) {
//       float tmp = t0;
//       t0 = t1;
//       t1 = tmp;
//     }
//     if (t0 > t_min)
//       t_min = t0;
//     if (t1 < t_max)
//       t_max = t1;
//     if (t_max < t_min)
//       return false;
//   }
//   return true;
// }

// AABB aabb_pad(const AABB *aabb) {
//   float delta = 1e-4f;
//   AABB padded;
//   for (int a = 0; a < 3; a++) {
//     if (aabb->x[a][1] - aabb->x[a][0] < delta) {
//       padded.x[a][0] = aabb->x[a][0] - delta;
//       padded.x[a][1] = aabb->x[a][1] + delta;
//     } else {
//       padded.x[a][0] = aabb->x[a][0];
//       padded.x[a][1] = aabb->x[a][1];
//     }
//   }
//   return padded;
// }

static bool Sphere_hit(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
static bool Quad_hit(const Quad *quad, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                          PCG32State *rng);
static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                        PCG32State *rng);
static bool ConstantMedium_hit(const ConstantMedium *constant_medium, const Ray *ray, float t_min, float t_max,
                               HitRecord *hit_record, PCG32State *rng);

bool Hittable_hit(Hittable obj, const Ray *ray, float t_min, float t_max, HitRecord *hit_record, PCG32State *rng) {
  switch (obj.type) {
  case HITTABLE_LIST:
    return HittableList_hit(obj.ptr, ray, t_min, t_max, hit_record, rng);
  case SPHERE:
    return Sphere_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case QUAD:
    return Quad_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case TRANSLATE:
    return Translate_hit(obj.ptr, ray, t_min, t_max, hit_record, rng);
  case ROTATE_Y:
    return RotateY_hit(obj.ptr, ray, t_min, t_max, hit_record, rng);
  case CONSTANT_MEDIUM:
    return ConstantMedium_hit(obj.ptr, ray, t_min, t_max, hit_record, rng);
  }
}

bool HittableList_hit(const HittableList *list, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                      PCG32State *rng) {
  bool hit_anything = false;

  for (int i = 0; i < list->size; i++)
    if (Hittable_hit(list->items[i], ray, t_min, t_max, hit_record, rng)) {
      t_max = hit_record->t;
      hit_anything = true;
    }

  return hit_anything;
}

static bool Sphere_hit(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  Vec3 oc = vec3_sub(ray->origin, sphere->center);
  float a = vec3_length2(ray->direction);
  float b = vec3_dot(oc, ray->direction);
  float c = vec3_length2(oc) - sphere->radius * sphere->radius;
  float disc = b * b - a * c;

  if (disc < 0)
    return false;

  float disc_sqrt = sqrtf(disc);
  float root = (-b - disc_sqrt) / a;
  if (root <= t_min || root >= t_max) {
    root = (-b + disc_sqrt) / a;
    if (root <= t_min || root >= t_max)
      return false;
  }

  hit_record->t = root;
  hit_record->p = ray_at(ray, root);

  Vec3 outward_normal = vec3_div(vec3_sub(hit_record->p, sphere->center), sphere->radius);
  hit_record->front_face = vec3_dot(ray->direction, outward_normal) < 0.0f;
  hit_record->normal = hit_record->front_face ? outward_normal : vec3_neg(outward_normal);
  hit_record->u = (atan2f(-outward_normal.x[2], outward_normal.x[0]) + (float)M_PI) * (float)M_1_PI * 0.5;
  hit_record->v = acosf(-outward_normal.x[1]) * M_1_PI;
  hit_record->material = sphere->material;

  return true;
}

static bool Quad_hit(const Quad *quad, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  float denom = vec3_dot(quad->normal, ray->direction);
  if (fabs(denom) < 1e-8f)
    return false;

  float t = (quad->D - vec3_dot(quad->normal, ray->origin)) / denom;
  if ((t < t_min) || (t > t_max))
    return false;

  Vec3 p = ray_at(ray, t);
  Vec3 planar_hitpt_vec = vec3_sub(p, quad->Q);
  float alpha = vec3_dot(quad->w, vec3_cross(planar_hitpt_vec, quad->v));
  float beta = vec3_dot(quad->w, vec3_cross(quad->u, planar_hitpt_vec));

  if ((alpha < 0) || (alpha > 1) || (beta < 0) || (beta > 1))
    return false;

  hit_record->u = alpha;
  hit_record->v = beta;
  hit_record->t = t;
  hit_record->p = p;
  hit_record->material = quad->material;
  hit_record->front_face = vec3_dot(ray->direction, quad->normal) < 0.0f;
  hit_record->normal = hit_record->front_face ? quad->normal : vec3_neg(quad->normal);

  return true;
}

static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                          PCG32State *rng) {
  Ray offset_r = {vec3_sub(ray->origin, translate->offset), ray->direction};

  if (!Hittable_hit(translate->object, &offset_r, t_min, t_max, hit_record, rng))
    return false;

  hit_record->p = vec3_add(hit_record->p, translate->offset);
  return true;
}

static Vec3 Vec3_rotate_y(Vec3 u, float cos_theta, float sin_theta) {
  return (Vec3){
      cos_theta * u.x[0] - sin_theta * u.x[2],
      u.x[1],
      sin_theta * u.x[0] + cos_theta * u.x[2],
  };
}

static Vec3 Vec3_rotate_y_inverse(Vec3 u, float cos_theta, float sin_theta) {
  return (Vec3){
      cos_theta * u.x[0] + sin_theta * u.x[2],
      u.x[1],
      -sin_theta * u.x[0] + cos_theta * u.x[2],
  };
}

static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *hit_record,
                        PCG32State *rng) {
  Ray rotated_r = {
      Vec3_rotate_y(ray->origin, rotate_y->cos_theta, rotate_y->sin_theta),
      Vec3_rotate_y(ray->direction, rotate_y->cos_theta, rotate_y->sin_theta),
  };

  if (!Hittable_hit(rotate_y->object, &rotated_r, t_min, t_max, hit_record, rng))
    return false;

  hit_record->p = Vec3_rotate_y_inverse(hit_record->p, rotate_y->cos_theta, rotate_y->sin_theta);
  hit_record->normal = Vec3_rotate_y_inverse(hit_record->normal, rotate_y->cos_theta, rotate_y->sin_theta);
  return true;
}

static bool ConstantMedium_hit(const ConstantMedium *constant_medium, const Ray *ray, float t_min, float t_max,
                               HitRecord *hit_record, PCG32State *rng) {
  HitRecord rec1, rec2;

  if (!Hittable_hit(constant_medium->boundary, ray, -INFINITY, INFINITY, &rec1, rng))
    return false;

  if (!Hittable_hit(constant_medium->boundary, ray, rec1.t + 0.0001f, INFINITY, &rec2, rng))
    return false;

  rec1.t = max(rec1.t, t_min);
  rec2.t = min(rec2.t, t_max);

  if (rec1.t >= rec2.t)
    return false;

  rec1.t = max(rec1.t, 0.0f);

  float ray_length = vec3_length(ray->direction);
  float distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
  float hit_distance = constant_medium->neg_inv_density * logf(pcg32_f32(rng));

  if (hit_distance > distance_inside_boundary)
    return false;

  // NOTE: we don't need to set normal and front_face, since Isotropic doesn't use them
  hit_record->t = rec1.t + hit_distance / ray_length;
  hit_record->p = ray_at(ray, hit_record->t);
  hit_record->material = constant_medium->phase_fn;
  return true;
}
