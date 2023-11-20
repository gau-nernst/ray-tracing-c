#include "hittable.h"
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
      {fminf(a.x[0][0], b.x[0][0]), fmaxf(a.x[0][1], b.x[0][1])},
      {fminf(a.x[1][0], b.x[1][0]), fmaxf(a.x[1][1], b.x[1][1])},
      {fminf(a.x[2][0], b.x[2][0]), fmaxf(a.x[2][1], b.x[2][1])},
  }};
}
static AABB AABB_pad(AABB bbox) {
  float delta = 1e-4f;
  AABB padded;
  for (int a = 0; a < 3; a++) {
    if (bbox.x[a][1] - bbox.x[a][0] < delta) {
      padded.x[a][0] = bbox.x[a][0] - delta;
      padded.x[a][1] = bbox.x[a][1] + delta;
    } else {
      padded.x[a][0] = bbox.x[a][0];
      padded.x[a][1] = bbox.x[a][1];
    }
  }
  return padded;
}

static AABB Hittable_bbox(Hittable obj) {
  switch (obj.type) {
  case HITTABLE_LIST:
    return ((HittableList *)obj.ptr)->bbox;
  case SPHERE:
    return ((Sphere *)obj.ptr)->bbox;
  case QUAD:
    return ((Quad *)obj.ptr)->bbox;
  case BVH_NODE:
    return ((BVHNode *)obj.ptr)->bbox;
  case TRANSLATE:
    return ((Translate *)obj.ptr)->bbox;
  case ROTATE_Y:
    return ((RotateY *)obj.ptr)->bbox;
  case CONSTANT_MEDIUM:
    return Hittable_bbox(((ConstantMedium *)obj.ptr)->boundary);
  default:
    assert(false && "Should not reached here");
  }
}

// for use with qsort()
typedef int (*Comparator)(const void *, const void *);
static int Hittable_compare_bbox(const Hittable *a, const Hittable *b, int axis) {
  float a_ = Hittable_bbox(*a).x[axis][0];
  float b_ = Hittable_bbox(*b).x[axis][0];
  return (a_ < b_) ? -1 : (a_ > b_) ? 1 : 0;
}
static int Hittable_compare_bbox_x(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 0); }
static int Hittable_compare_bbox_y(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 1); }
static int Hittable_compare_bbox_z(const void *a, const void *b) { return Hittable_compare_bbox(a, b, 2); }

void HittableList_init(HittableList *list, size_t max_size) {
  *list = (HittableList){max_size, 0, my_malloc(sizeof(list->items[0]) * max_size), AABB_EMPTY};
}
HittableList *HittableList_new(size_t max_size) define_init_new(HittableList, max_size);
void HittableList_append(HittableList *list, Hittable item) {
  assert((list->size < list->max_size) && "List is full");
  list->items[list->size++] = item;
  list->bbox = AABB_from_AABB(list->bbox, Hittable_bbox(item));
}

void Sphere_init(Sphere *sphere, Vec3 center, float radius, Material mat) {
  *sphere = (Sphere){center, radius, mat, AABB_from_Vec3(vec3_sub(center, radius), vec3_add(center, radius))};
}
Sphere *Sphere_new(Vec3 center, float radius, Material mat) define_init_new(Sphere, center, radius, mat);

void Quad_init(Quad *quad, Vec3 Q, Vec3 u, Vec3 v, Material mat) {
  quad->Q = Q;
  quad->u = u;
  quad->v = v;
  quad->material = mat;
  quad->bbox = AABB_pad(AABB_from_Vec3(Q, vec3_add(Q, u, v)));

  Vec3 n = vec3_cross(quad->u, quad->v);
  quad->normal = vec3_normalize(n);
  quad->D = vec3_dot(quad->normal, quad->Q);
  quad->w = vec3_div(n, vec3_length2(n));
}
Quad *Quad_new(Vec3 Q, Vec3 u, Vec3 v, Material mat) define_init_new(Quad, Q, u, v, mat);

HittableList *Box_new(Vec3 a, Vec3 b, Material mat) {
  HittableList *list = HittableList_new(6);

  Vec3 min_p = vec3_min(a, b);
  Vec3 max_p = vec3_max(a, b);

  Vec3 dx = vec3(max_p.x - min_p.x, 0, 0);
  Vec3 dy = vec3(0, max_p.y - min_p.y, 0);
  Vec3 dz = vec3(0, 0, max_p.z - min_p.z);

  HittableList_append(list, hittable(Quad_new(vec3(min_p.x, min_p.y, max_p.z), dx, dy, mat)));           // front
  HittableList_append(list, hittable(Quad_new(vec3(max_p.x, min_p.y, max_p.z), vec3_neg(dz), dy, mat))); // right
  HittableList_append(list, hittable(Quad_new(vec3(max_p.x, min_p.y, min_p.z), vec3_neg(dx), dy, mat))); // back
  HittableList_append(list, hittable(Quad_new(vec3(min_p.x, min_p.y, min_p.z), dz, dy, mat)));           // left
  HittableList_append(list, hittable(Quad_new(vec3(min_p.x, max_p.y, max_p.z), dx, vec3_neg(dz), mat))); // top
  HittableList_append(list, hittable(Quad_new(vec3(min_p.x, min_p.y, min_p.z), dx, dz, mat)));           // bottom

  return list;
}

static void _BVHNode_init(BVHNode *bvh, const Hittable *list, size_t n, PCG32State *rng) {
  // make a copy
  Hittable *list_ = my_malloc(sizeof(Hittable) * n);
  for (int i = 0; i < n; i++)
    list_[i] = list[i];

  int axis = pcg32_u32_between(rng, 0, 3);

  if (n == 1) {
    bvh->left = list_[0];
    bvh->right = list_[0];
  } else if (n == 2) {
    if (Hittable_compare_bbox(list_, list_ + 1, axis) < 0) {
      bvh->left = list_[0];
      bvh->right = list_[1];
    } else {
      bvh->left = list_[1];
      bvh->right = list_[0];
    }
  } else {
    Comparator comparator = (axis == 0)   ? Hittable_compare_bbox_x
                            : (axis == 1) ? Hittable_compare_bbox_y
                                          : Hittable_compare_bbox_z;
    qsort(list_, n, sizeof(Hittable), comparator);

    size_t mid = n / 2;
    // bvh->left = hittable(BVHNode_new(list_, mid, rng));
    // bvh->right = hittable(BVHNode_new(list_ + mid, n - mid, rng));

    BVHNode *left = my_malloc(sizeof(BVHNode));
    _BVHNode_init(left, list_, mid, rng);
    bvh->left = hittable(left);

    BVHNode *right = my_malloc(sizeof(BVHNode));
    _BVHNode_init(right, list_ + mid, n - mid, rng);
    bvh->right = hittable(right);
  }
  free(list_); // we don't need this anymore
  bvh->bbox = AABB_from_AABB(Hittable_bbox(bvh->left), Hittable_bbox(bvh->right));
}
void BVHNode_init(BVHNode *bvh, const HittableList *list, PCG32State *rng) {
  _BVHNode_init(bvh, list->items, list->size, rng);
}
BVHNode *BVHNode_new(const HittableList *list, PCG32State *rng) define_init_new(BVHNode, list, rng);

void Translate_init(Translate *translate, Hittable object, Vec3 offset) {
  AABB bbox = Hittable_bbox(object);
  for (int i = 0; i < 3; i++) {
    bbox.x[i][0] += offset.values[i];
    bbox.x[i][1] += offset.values[i];
  }
  *translate = (Translate){object, offset, bbox};
}
Translate *Translate_new(Hittable object, Vec3 offset) define_init_new(Translate, object, offset);

void RotateY_init(RotateY *rotate_y, Hittable object, float angle) {
  angle = angle * (float)M_PI / 180.f;
  float sin_theta = sinf(angle);
  float cos_theta = cosf(angle);

  AABB bbox = Hittable_bbox(object);
  Vec3 min_p = vec3(INFINITY, INFINITY, INFINITY);
  Vec3 max_p = vec3(-INFINITY, -INFINITY, -INFINITY);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        Vec3 tester = vec3((float)i * bbox.x[0][1] + (float)(1 - i) * bbox.x[0][0],
                           (float)j * bbox.x[1][1] + (float)(1 - j) * bbox.x[1][0],
                           (float)k * bbox.x[2][1] + (float)(1 - k) * bbox.x[2][0]);
        tester = Vec3_rotate_y_inverse(tester, cos_theta, sin_theta);
        min_p = vec3_min(min_p, tester);
        max_p = vec3_max(max_p, tester);
      }

  *rotate_y = (RotateY){object, sin_theta, cos_theta, AABB_from_Vec3(min_p, max_p)};
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

static bool AABB_hit(const AABB *aabb, const Ray *ray, float t_min, float t_max) {
  for (int i = 0; i < 3; i++) {
    float invD = 1.0f / ray->direction.values[i];
    float t0 = (aabb->x[i][0] - ray->origin.values[i]) * invD;
    float t1 = (aabb->x[i][1] - ray->origin.values[i]) * invD;

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

static bool Sphere_hit(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *rec);
static bool Quad_hit(const Quad *quad, const Ray *ray, float t_min, float t_max, HitRecord *rec);
static bool BVHNode_hit(const BVHNode *bvh, const Ray *ray, float t_min, float t_max, HitRecord *rec, PCG32State *rng);
static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                          PCG32State *rng);
static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                        PCG32State *rng);
static bool ConstantMedium_hit(const ConstantMedium *constant_medium, const Ray *ray, float t_min, float t_max,
                               HitRecord *rec, PCG32State *rng);

bool Hittable_hit(Hittable obj, const Ray *ray, float t_min, float t_max, HitRecord *rec, PCG32State *rng) {
  switch (obj.type) {
  case HITTABLE_LIST:
    return HittableList_hit(obj.ptr, ray, t_min, t_max, rec, rng);
  case SPHERE:
    return Sphere_hit(obj.ptr, ray, t_min, t_max, rec);
  case QUAD:
    return Quad_hit(obj.ptr, ray, t_min, t_max, rec);
  case BVH_NODE:
    return BVHNode_hit(obj.ptr, ray, t_min, t_max, rec, rng);
  case TRANSLATE:
    return Translate_hit(obj.ptr, ray, t_min, t_max, rec, rng);
  case ROTATE_Y:
    return RotateY_hit(obj.ptr, ray, t_min, t_max, rec, rng);
  case CONSTANT_MEDIUM:
    return ConstantMedium_hit(obj.ptr, ray, t_min, t_max, rec, rng);
  }
}

bool HittableList_hit(const HittableList *list, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                      PCG32State *rng) {
  bool hit_anything = false;

  for (int i = 0; i < list->size; i++)
    if (Hittable_hit(list->items[i], ray, t_min, t_max, rec, rng)) {
      t_max = rec->t;
      hit_anything = true;
    }

  return hit_anything;
}

static bool Sphere_hit(const Sphere *sphere, const Ray *ray, float t_min, float t_max, HitRecord *rec) {
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

  rec->t = root;
  rec->p = ray_at(ray, root);

  Vec3 outward_normal = vec3_div(vec3_sub(rec->p, sphere->center), sphere->radius);
  rec->front_face = vec3_dot(ray->direction, outward_normal) < 0.0f;
  rec->normal = rec->front_face ? outward_normal : vec3_neg(outward_normal);
  rec->u = (atan2f(-outward_normal.z, outward_normal.x) + (float)M_PI) * (float)M_1_PI * 0.5;
  rec->v = acosf(-outward_normal.y) * M_1_PI;
  rec->material = sphere->material;

  return true;
}

static bool Quad_hit(const Quad *quad, const Ray *ray, float t_min, float t_max, HitRecord *rec) {
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

  rec->u = alpha;
  rec->v = beta;
  rec->t = t;
  rec->p = p;
  rec->material = quad->material;
  rec->front_face = vec3_dot(ray->direction, quad->normal) < 0.0f;
  rec->normal = rec->front_face ? quad->normal : vec3_neg(quad->normal);

  return true;
}

static bool BVHNode_hit(const BVHNode *bvh, const Ray *ray, float t_min, float t_max, HitRecord *rec, PCG32State *rng) {
  if (!AABB_hit(&bvh->bbox, ray, t_min, t_max))
    return false;

  // NOTE: we need to check for both left and right, since we don't know which one is closer.
  bool hit_left = Hittable_hit(bvh->left, ray, t_min, t_max, rec, rng);
  if (hit_left)
    t_max = rec->t;
  bool hit_right = Hittable_hit(bvh->right, ray, t_min, t_max, rec, rng);
  return hit_left || hit_right;
}

static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                          PCG32State *rng) {
  Ray offset_r = {vec3_sub(ray->origin, translate->offset), ray->direction};

  if (!Hittable_hit(translate->object, &offset_r, t_min, t_max, rec, rng))
    return false;

  rec->p = vec3_add(rec->p, translate->offset);
  return true;
}

static Vec3 Vec3_rotate_y(Vec3 u, float cos_theta, float sin_theta) {
  return vec3(cos_theta * u.x - sin_theta * u.z, u.y, sin_theta * u.x + cos_theta * u.z);
}

static Vec3 Vec3_rotate_y_inverse(Vec3 u, float cos_theta, float sin_theta) {
  return vec3(cos_theta * u.x + sin_theta * u.z, u.y, -sin_theta * u.x + cos_theta * u.z);
}

static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *rec,
                        PCG32State *rng) {
  Ray rotated_r = {
      Vec3_rotate_y(ray->origin, rotate_y->cos_theta, rotate_y->sin_theta),
      Vec3_rotate_y(ray->direction, rotate_y->cos_theta, rotate_y->sin_theta),
  };

  if (!Hittable_hit(rotate_y->object, &rotated_r, t_min, t_max, rec, rng))
    return false;

  rec->p = Vec3_rotate_y_inverse(rec->p, rotate_y->cos_theta, rotate_y->sin_theta);
  rec->normal = Vec3_rotate_y_inverse(rec->normal, rotate_y->cos_theta, rotate_y->sin_theta);
  return true;
}

static bool ConstantMedium_hit(const ConstantMedium *constant_medium, const Ray *ray, float t_min, float t_max,
                               HitRecord *rec, PCG32State *rng) {
  HitRecord rec1, rec2;

  if (!Hittable_hit(constant_medium->boundary, ray, -INFINITY, INFINITY, &rec1, rng))
    return false;

  if (!Hittable_hit(constant_medium->boundary, ray, rec1.t + 0.0001f, INFINITY, &rec2, rng))
    return false;

  rec1.t = fmaxf(rec1.t, t_min);
  rec2.t = fminf(rec2.t, t_max);

  if (rec1.t >= rec2.t)
    return false;

  rec1.t = max(rec1.t, 0.0f);

  float ray_length = vec3_length(ray->direction);
  float distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
  float hit_distance = constant_medium->neg_inv_density * logf(pcg32_f32(rng));

  if (hit_distance > distance_inside_boundary)
    return false;

  // NOTE: we don't need to set normal and front_face, since Isotropic doesn't use them
  rec->t = rec1.t + hit_distance / ray_length;
  rec->p = ray_at(ray, rec->t);
  rec->material = constant_medium->phase_fn;
  return true;
}
