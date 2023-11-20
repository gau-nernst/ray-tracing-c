#include "hittable.h"
#include "pcg32.h"
#include "vec3.h"
#include <math.h>

const static AABB AABB_EMPTY = {{{INFINITY, -INFINITY}, {INFINITY, -INFINITY}, {INFINITY, -INFINITY}}};
static AABB AABB_from_Vec3(Vec3 a, Vec3 b);
static AABB AABB_from_AABB(AABB a, AABB b);
static AABB AABB_pad(AABB bbox);
static Vec3 Vec3_rotate_y(Vec3 u, float cos_theta, float sin_theta);
static Vec3 Vec3_rotate_y_inverse(Vec3 u, float cos_theta, float sin_theta);

Vec3 ray_at(const Ray *ray, float t) { return vec3vec3_add(ray->origin, vec3float_mul(ray->direction, t)); }

static AABB AABB_from_Vec3(Vec3 a, Vec3 b) {
  return (AABB){{
      {fminf(a.x, b.x), fmaxf(a.x, b.x)},
      {fminf(a.y, b.y), fmaxf(a.y, b.y)},
      {fminf(a.z, b.z), fmaxf(a.z, b.z)},
  }};
}
static AABB AABB_from_AABB(AABB a, AABB b) {
  return (AABB){{
      {fminf(a.x[0], b.x[0]), fmaxf(a.x[1], b.x[1])},
      {fminf(a.y[0], b.y[0]), fmaxf(a.y[1], b.y[1])},
      {fminf(a.z[0], b.z[0]), fmaxf(a.z[1], b.z[1])},
  }};
}
static AABB AABB_pad(AABB bbox) {
  float delta = 1e-4f;
  AABB padded;
  for (int a = 0; a < 3; a++) {
    if (bbox.values[a][1] - bbox.values[a][0] < delta) {
      padded.values[a][0] = bbox.values[a][0] - delta;
      padded.values[a][1] = bbox.values[a][1] + delta;
    } else {
      padded.values[a][0] = bbox.values[a][0];
      padded.values[a][1] = bbox.values[a][1];
    }
  }
  return padded;
}

// static AABB Hittable_bbox(Hittable obj) {
//   switch (obj.type) {
//   case HITTABLE_LIST:
//     return ((HittableList *)obj.ptr)->bbox;
//   case SPHERE:
//     return ((Sphere *)obj.ptr)->bbox;
//   case QUAD:
//     return ((Quad *)obj.ptr)->bbox;
//   case BVH_NODE:
//     return ((BVHNode *)obj.ptr)->bbox;
//   case TRANSLATE:
//     return ((Translate *)obj.ptr)->bbox;
//   case ROTATE_Y:
//     return ((RotateY *)obj.ptr)->bbox;
//   case CONSTANT_MEDIUM:
//     return Hittable_bbox(((ConstantMedium *)obj.ptr)->boundary);
//   default:
//     assert(false && "Should not reached here");
//   }
// }

// for use with qsort()
typedef int (*Comparator)(const void *, const void *);
static int Hittable_compare_bbox(const Hittable **a, const Hittable **b, int axis) {
  float a_ = (*a)->bbox(*a).values[axis][0];
  float b_ = (*b)->bbox(*b).values[axis][0];
  return (a_ < b_) ? -1 : (a_ > b_) ? 1 : 0;
}
static int Hittable_compare_bbox_x(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 0); }
static int Hittable_compare_bbox_y(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 1); }
static int Hittable_compare_bbox_z(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 2); }

static HittableHitFn HittableList_hit;
static AABB HittableList_bbox(const Hittable *self_) { return ((HittableList *)self_)->bbox; }
void HittableList_init(HittableList *self, size_t max_size) {
  *self = (HittableList){
      {HittableList_hit, HittableList_bbox}, max_size, 0, my_malloc(sizeof(Hittable *) * max_size), AABB_EMPTY};
}
Hittable *HittableList_new(size_t max_size) define_init_new(HittableList, max_size);
void HittableList_append(HittableList *self, Hittable *item) {
  assert((self->size < self->max_size) && "List is full");
  self->items[self->size++] = item;
  self->bbox = AABB_from_AABB(self->bbox, item->bbox(item));
}
static bool HittableList_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                             PCG32State *rng) {
  HittableList *self = (HittableList *)self_;
  bool hit_anything = false;

  for (int i = 0; i < self->size; i++) {
    Hittable *item = self->items[i];
    if (item->hit(item, ray, t_min, t_max, rec, rng)) {
      t_max = rec->t;
      hit_anything = true;
    }
  }

  return hit_anything;
}

static HittableHitFn Sphere_hit;
static AABB Sphere_bbox(const Hittable *self_) { return ((Sphere *)self_)->bbox; }
void Sphere_init(Sphere *sphere, Vec3 center, float radius, Material mat) {
  *sphere = (Sphere){{Sphere_hit, Sphere_bbox},
                     center,
                     radius,
                     mat,
                     AABB_from_Vec3(vec3_sub(center, radius), vec3_add(center, radius))};
}
Hittable *Sphere_new(Vec3 center, float radius, Material mat) define_init_new(Sphere, center, radius, mat);

static bool Sphere_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                       PCG32State *rng) {
  Sphere *self = (Sphere *)self_;

  Vec3 oc = vec3_sub(ray->origin, self->center);
  float a = vec3_length2(ray->direction);
  float b = vec3_dot(oc, ray->direction);
  float c = vec3_length2(oc) - self->radius * self->radius;
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

  rec->t = root;
  rec->p = ray_at(ray, root);

  Vec3 outward_normal = vec3_div(vec3_sub(rec->p, self->center), self->radius);
  rec->front_face = vec3_dot(ray->direction, outward_normal) < 0.0f;
  rec->normal = rec->front_face ? outward_normal : vec3_neg(outward_normal);
  rec->u = (atan2f(-outward_normal.z, outward_normal.x) + (float)M_PI) * (float)M_1_PI * 0.5;
  rec->v = acosf(-outward_normal.y) * M_1_PI;
  rec->material = self->material;

  return true;
}

static HittableHitFn Quad_hit;
static AABB Quad_bbox(const Hittable *self_) { return ((Quad *)self_)->bbox; }
void Quad_init(Quad *self, Vec3 Q, Vec3 u, Vec3 v, Material mat) {
  self->hittable = (Hittable){Quad_hit, Quad_bbox};
  self->Q = Q;
  self->u = u;
  self->v = v;
  self->material = mat;
  self->bbox = AABB_pad(AABB_from_Vec3(Q, vec3_add(Q, u, v)));

  Vec3 n = vec3_cross(self->u, self->v);
  self->normal = vec3_normalize(n);
  self->D = vec3_dot(self->normal, self->Q);
  self->w = vec3_div(n, vec3_length2(n));
}
Hittable *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material mat) define_init_new(Quad, Q, u, v, mat);
static bool Quad_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec, PCG32State *rng) {
  Quad *self = (Quad *)self_;

  float denom = vec3_dot(self->normal, ray->direction);
  if (fabs(denom) < 1e-8f)
    return false;

  float t = (self->D - vec3_dot(self->normal, ray->origin)) / denom;
  if ((t < t_min) || (t > t_max))
    return false;

  Vec3 p = ray_at(ray, t);
  Vec3 planar_hitpt_vec = vec3_sub(p, self->Q);
  float alpha = vec3_dot(self->w, vec3_cross(planar_hitpt_vec, self->v));
  float beta = vec3_dot(self->w, vec3_cross(self->u, planar_hitpt_vec));

  if ((alpha < 0) || (alpha > 1) || (beta < 0) || (beta > 1))
    return false;

  rec->u = alpha;
  rec->v = beta;
  rec->t = t;
  rec->p = p;
  rec->material = self->material;
  rec->front_face = vec3_dot(ray->direction, self->normal) < 0.0f;
  rec->normal = rec->front_face ? self->normal : vec3_neg(self->normal);

  return true;
}

Hittable *Box_new(Vec3 a, Vec3 b, Material mat) {
  HittableList *list = (HittableList *)HittableList_new(6);

  Vec3 min_p = vec3_min(a, b);
  Vec3 max_p = vec3_max(a, b);

  Vec3 dx = vec3(max_p.x - min_p.x, 0, 0);
  Vec3 dy = vec3(0, max_p.y - min_p.y, 0);
  Vec3 dz = vec3(0, 0, max_p.z - min_p.z);

  HittableList_append(list, Quad_new(vec3(min_p.x, min_p.y, max_p.z), dx, dy, mat));           // front
  HittableList_append(list, Quad_new(vec3(max_p.x, min_p.y, max_p.z), vec3_neg(dz), dy, mat)); // right
  HittableList_append(list, Quad_new(vec3(max_p.x, min_p.y, min_p.z), vec3_neg(dx), dy, mat)); // back
  HittableList_append(list, Quad_new(vec3(min_p.x, min_p.y, min_p.z), dz, dy, mat));           // left
  HittableList_append(list, Quad_new(vec3(min_p.x, max_p.y, max_p.z), dx, vec3_neg(dz), mat)); // top
  HittableList_append(list, Quad_new(vec3(min_p.x, min_p.y, min_p.z), dx, dz, mat));           // bottom

  return (Hittable *)list;
}

static HittableHitFn BVHNode_hit;
static AABB BVHNode_bbox(const Hittable *self_) { return ((BVHNode *)self_)->bbox; }
static void _BVHNode_init(BVHNode *self, Hittable **list_, size_t n, PCG32State *rng);
void BVHNode_init(BVHNode *self, const HittableList *list, PCG32State *rng) {
  _BVHNode_init(self, list->items, list->size, rng);
}
Hittable *BVHNode_new(const HittableList *list, PCG32State *rng) define_init_new(BVHNode, list, rng);
static void _BVHNode_init(BVHNode *self, Hittable **list_, size_t n, PCG32State *rng) {
  self->hittable = (Hittable){BVHNode_hit, BVHNode_bbox};

  // make a copy
  Hittable **list = my_malloc(sizeof(Hittable *) * n);
  for (int i = 0; i < n; i++)
    list[i] = list_[i];

  int axis = pcg32_u32_between(rng, 0, 3);

  if (n == 1) {
    self->left = list[0];
    self->right = list[0];
  } else if (n == 2) {
    if (Hittable_compare_bbox(list, list + 1, axis) < 0) {
      self->left = list[0];
      self->right = list[1];
    } else {
      self->left = list[1];
      self->right = list[0];
    }
  } else {
    Comparator comparator = (axis == 0)   ? Hittable_compare_bbox_x
                            : (axis == 1) ? Hittable_compare_bbox_y
                                          : Hittable_compare_bbox_z;
    qsort(list, n, sizeof(Hittable), comparator);

    size_t mid = n / 2;

    BVHNode *left = my_malloc(sizeof(BVHNode));
    _BVHNode_init(left, list, mid, rng);
    self->left = (Hittable *)left;

    BVHNode *right = my_malloc(sizeof(BVHNode));
    _BVHNode_init(right, list + mid, n - mid, rng);
    self->right = (Hittable *)right;
  }
  free(list); // we don't need this anymore
  self->bbox = AABB_from_AABB(self->left->bbox(self->left), self->right->bbox(self->right));
}
static bool AABB_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max) {
  for (int i = 0; i < 3; i++) {
    float invD = 1.0f / ray->direction.values[i];
    float t0 = (aabb->values[i][0] - ray->origin.values[i]) * invD;
    float t1 = (aabb->values[i][1] - ray->origin.values[i]) * invD;

    if (invD < 0) {
      float tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    t_min = fmaxf(t_min, t0);
    t_max = fminf(t_max, t1);
    if (t_max <= t_min)
      return false;
  }
  return true;
}
static bool BVHNode_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                        PCG32State *rng) {
  BVHNode *self = (BVHNode *)self_;
  if (!AABB_hit(&self->bbox, ray, t_min, t_max))
    return false;

  // NOTE: we need to check for both left and right, since we don't know which one is closer.
  bool hit_left = self->left->hit(self->left, ray, t_min, t_max, rec, rng);
  if (hit_left)
    t_max = rec->t;
  bool hit_right = self->right->hit(self->right, ray, t_min, t_max, rec, rng);
  return hit_left || hit_right;
}

static HittableHitFn Translate_hit;
static AABB Translate_bbox(const Hittable *self) { return ((Translate *)self)->bbox; }
void Translate_init(Translate *self, Hittable *object, Vec3 offset) {
  AABB bbox = object->bbox(object);
  for (int i = 0; i < 2; i++) {
    bbox.x[i] += offset.x;
    bbox.y[i] += offset.y;
    bbox.z[i] += offset.z;
  }
  *self = (Translate){{Translate_hit, Translate_bbox}, object, offset, bbox};
}
Hittable *Translate_new(Hittable *object, Vec3 offset) define_init_new(Translate, object, offset);
static bool Translate_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                          PCG32State *rng) {
  Translate *self = (Translate *)self_;
  Ray offset_r = {vec3_sub(ray->origin, self->offset), ray->direction};

  if (!self->object->hit(self->object, &offset_r, t_min, t_max, rec, rng))
    return false;

  rec->p = vec3_add(rec->p, self->offset);
  return true;
}

static HittableHitFn RotateY_hit;
static AABB RotateY_bbox(const Hittable *self) { return ((RotateY *)self)->bbox; }
void RotateY_init(RotateY *self, Hittable *object, float angle) {
  angle = angle * (float)M_PI / 180.f;
  float sin_theta = sinf(angle);
  float cos_theta = cosf(angle);

  AABB bbox = object->bbox(object);
  Vec3 min_p = vec3(INFINITY, INFINITY, INFINITY);
  Vec3 max_p = vec3(-INFINITY, -INFINITY, -INFINITY);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        Vec3 tester =
            vec3((float)i * bbox.x[1] + (float)(1 - i) * bbox.x[0], (float)j * bbox.y[1] + (float)(1 - j) * bbox.y[0],
                 (float)k * bbox.z[1] + (float)(1 - k) * bbox.z[0]);
        tester = Vec3_rotate_y_inverse(tester, cos_theta, sin_theta);
        min_p = vec3_min(min_p, tester);
        max_p = vec3_max(max_p, tester);
      }

  *self = (RotateY){{RotateY_hit, RotateY_bbox}, object, sin_theta, cos_theta, AABB_from_Vec3(min_p, max_p)};
}
Hittable *RotateY_new(Hittable *object, float angle) define_init_new(RotateY, object, angle);
static bool RotateY_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                        PCG32State *rng) {
  RotateY *self = (RotateY *)self_;
  Ray rotated_r = {
      Vec3_rotate_y(ray->origin, self->cos_theta, self->sin_theta),
      Vec3_rotate_y(ray->direction, self->cos_theta, self->sin_theta),
  };

  if (!self->object->hit(self->object, &rotated_r, t_min, t_max, rec, rng))
    return false;

  rec->p = Vec3_rotate_y_inverse(rec->p, self->cos_theta, self->sin_theta);
  rec->normal = Vec3_rotate_y_inverse(rec->normal, self->cos_theta, self->sin_theta);
  return true;
}

static HittableHitFn ConstantMedium_hit;
static AABB ConstantMedium_bbox(const Hittable *self) {
  Hittable *boundary = ((ConstantMedium *)self)->boundary;
  return boundary->bbox(boundary);
}
void ConstantMedium_init(ConstantMedium *self, Hittable *boundary, float density, Texture *albedo) {
  *self = (ConstantMedium){
      {ConstantMedium_hit, ConstantMedium_bbox},
      boundary,
      -1.0f / density,
      material(Isotropic_new(albedo)),
  };
}
Hittable *ConstantMedium_new(Hittable *boundary, float density, Texture *albedo)
    define_init_new(ConstantMedium, boundary, density, albedo);
static bool ConstantMedium_hit(const Hittable *self_, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                               PCG32State *rng) {
  ConstantMedium *self = (ConstantMedium *)self_;
  HitRecord rec1, rec2;

  if (!self->boundary->hit(self->boundary, ray, -INFINITY, INFINITY, &rec1, rng))
    return false;

  if (!self->boundary->hit(self->boundary, ray, rec1.t + 0.0001f, INFINITY, &rec2, rng))
    return false;

  rec1.t = fmaxf(rec1.t, t_min);
  rec2.t = fminf(rec2.t, t_max);

  if (rec1.t >= rec2.t)
    return false;

  rec1.t = max(rec1.t, 0.0f);

  float ray_length = vec3_length(ray->direction);
  float distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
  float hit_distance = self->neg_inv_density * logf(pcg32_f32(rng));

  if (hit_distance > distance_inside_boundary)
    return false;

  // NOTE: we don't need to set normal and front_face, since Isotropic doesn't use them
  rec->t = rec1.t + hit_distance / ray_length;
  rec->p = ray_at(ray, rec->t);
  rec->material = self->phase_fn;
  return true;
}

static Vec3 Vec3_rotate_y(Vec3 u, float cos_theta, float sin_theta) {
  return vec3(cos_theta * u.x - sin_theta * u.z, u.y, sin_theta * u.x + cos_theta * u.z);
}

static Vec3 Vec3_rotate_y_inverse(Vec3 u, float cos_theta, float sin_theta) {
  return vec3(cos_theta * u.x + sin_theta * u.z, u.y, -sin_theta * u.x + cos_theta * u.z);
}
