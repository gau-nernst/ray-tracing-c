#include "raytracing.h"
#include "utils.h"
#include <stdio.h>

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
static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);
static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *hit_record);

bool Hittable_hit(Hittable obj, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  switch (obj.type) {
  case HITTABLE_LIST:
    return HittableList_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case SPHERE:
    return Sphere_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case QUAD:
    return Quad_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case TRANSLATE:
    return Translate_hit(obj.ptr, ray, t_min, t_max, hit_record);
  case ROTATE_Y:
    return RotateY_hit(obj.ptr, ray, t_min, t_max, hit_record);
  }
}

bool HittableList_hit(const HittableList *list, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  bool hit_anything = false;

  for (int i = 0; i < list->size; i++)
    if (Hittable_hit(list->items[i], ray, t_min, t_max, hit_record)) {
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

static bool Translate_hit(const Translate *translate, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  Ray offset_r = {vec3_sub(ray->origin, translate->offset), ray->direction};

  if (!Hittable_hit(translate->object, &offset_r, t_min, t_max, hit_record))
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

static bool RotateY_hit(const RotateY *rotate_y, const Ray *ray, float t_min, float t_max, HitRecord *hit_record) {
  Ray rotated_r = {
      Vec3_rotate_y(ray->origin, rotate_y->cos_theta, rotate_y->sin_theta),
      Vec3_rotate_y(ray->direction, rotate_y->cos_theta, rotate_y->sin_theta),
  };

  if (!Hittable_hit(rotate_y->object, &rotated_r, t_min, t_max, hit_record))
    return false;

  hit_record->p = Vec3_rotate_y_inverse(hit_record->p, rotate_y->cos_theta, rotate_y->sin_theta);
  hit_record->normal = Vec3_rotate_y_inverse(hit_record->normal, rotate_y->cos_theta, rotate_y->sin_theta);
  return true;
}

void World_init(World *world, size_t max_objects, size_t max_materials) {
  HittableList_init(&world->objects, max_objects);
  MaterialList_init(&world->materials, max_materials);
}

void Camera_init(Camera *camera) {
  camera->img_height = (int)((float)camera->img_width / camera->aspect_ratio);

  float viewport_height = 2.0 * tanf(camera->vfov * (float)M_PI / 360.0f) * camera->focal_length;
  float viewport_width = viewport_height * (float)camera->img_width / (float)camera->img_height;

  camera->w = vec3_unit(vec3_sub(camera->look_from, camera->look_to));
  camera->u = vec3_cross(camera->vup, camera->w);
  camera->v = vec3_cross(camera->w, camera->u);

  Vec3 viewport_u = vec3_mul(camera->u, viewport_width);   // scan from left to right
  Vec3 viewport_v = vec3_mul(camera->v, -viewport_height); // scan from top to bottom

  camera->pixel_delta_u = vec3_div(viewport_u, (float)camera->img_width);
  camera->pixel_delta_v = vec3_div(viewport_v, (float)camera->img_height);

  Vec3 viewport_upper_left = vec3_add(camera->look_from, vec3_mul(camera->w, -camera->focal_length),
                                      vec3_mul(viewport_u, -0.5f), vec3_mul(viewport_v, -0.5f));
  camera->pixel00_loc =
      vec3_add(viewport_upper_left, vec3_mul(camera->pixel_delta_u, 0.5f), vec3_mul(camera->pixel_delta_v, 0.5f));

  float dof_radius = camera->focal_length * tanf(camera->dof_angle * (float)M_PI / 360.0f);
  camera->dof_disc_u = vec3_mul(camera->u, dof_radius);
  camera->dof_disc_v = vec3_mul(camera->v, dof_radius);
}

static Vec3 Camera_ray_color(const Camera *camera, const Ray *ray, const World *world, int depth, PCG32State *rng) {
  if (depth <= 0)
    return (Vec3){0, 0, 0};

  HitRecord hit_record;
  if (HittableList_hit(&world->objects, ray, 1e-3f, INFINITY, &hit_record)) {
    Ray new_ray;
    new_ray.origin = hit_record.p;
    Vec3 scatter_color;

    if (scatter(ray->direction, &hit_record, rng, &new_ray.direction, &scatter_color))
      scatter_color =
          vec3_mul(Camera_ray_color(camera, &new_ray, world, depth - 1, rng), scatter_color); // spawn new ray

    return vec3_add(scatter_color, emit(&hit_record));
  }

  // scene background
  return camera->background;

  // old background
  // Vec3 direction = vec3_unit(ray->direction);
  // float a = 0.5f * (direction.y + 1.0f); // [-1,1] -> [0,1]

  // Vec3 WHITE = {1.0f, 1.0f, 1.0f};
  // Vec3 BLUE = {0.5f, 0.7f, 1.0f};
  // return vec3_lerp(WHITE, BLUE, a);
}

void Camera_render(const Camera *camera, const World *world, uint8_t *buffer) {
  for (int j = 0; j < camera->img_height; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", camera->img_height - j);

    int i; // C89 for loop for MSVC OpenMP
#pragma omp parallel for private(i) schedule(static, 1)
    for (i = 0; i < camera->img_width; i++) {
      PCG32State rng;
      pcg32_seed(&rng, 17 + j, 23 + i);

      Vec3 pixel_pos = vec3_add(camera->pixel00_loc, vec3_mul(camera->pixel_delta_u, (float)i),
                                vec3_mul(camera->pixel_delta_v, (float)j));
      Vec3 pixel_color = {0, 0, 0};

      for (int sample = 0; sample < camera->samples_per_pixel; sample++) {
        // square sampling
        // another option: sinc sampling
        // TODO: use 64-bit PRNG to generate 2 numbers at once
        float px = pcg32_f32_between(&rng, -0.5f, 0.5f);
        float py = pcg32_f32_between(&rng, -0.5f, 0.5f);

        Ray ray;
        if (camera->dof_angle > 0.0f) {
          // sample points around camera.look_from like a (thin) lens/aperture
          float a, b;
          for (;;) {
            a = pcg32_f32_between(&rng, -1.0f, 1.0f);
            b = pcg32_f32_between(&rng, -1.0f, 1.0f);
            if (a * a + b * b < 1.0f)
              break;
          }
          ray.origin = vec3_add(camera->look_from, vec3_mul(camera->dof_disc_u, a), vec3_mul(camera->dof_disc_v, b));
        } else {
          ray.origin = camera->look_from;
        }
        ray.direction = vec3_add(pixel_pos, vec3_mul(camera->pixel_delta_u, px), vec3_mul(camera->pixel_delta_v, py),
                                 vec3_neg(ray.origin));

        pixel_color = vec3_add(pixel_color, Camera_ray_color(camera, &ray, world, camera->max_depth, &rng));
      }

      pixel_color = vec3_div(pixel_color, (float)camera->samples_per_pixel);
      for (int c = 0; c < 3; c++)
        buffer[(j * camera->img_width + i) * 3 + c] = (int)(256.0f * clamp(sqrtf(pixel_color.x[c]), 0.0f, 0.999f));
    }
  }
  fprintf(stderr, "\nDone\n");
}
