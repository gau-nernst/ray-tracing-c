#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdlib.h>

// gcc with -std=c11 doesn't define M_PI
#define M_PI 3.14159265358979323846264338327950288
#define M_1_PI 0.318309886183790671537767526745028724

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#define clamp(x, lo, hi) min(max(x, lo), hi)

#define define_struct_new(Type, ...)                                                                                   \
  {                                                                                                                    \
    Type *obj = my_malloc(sizeof(Type));                                                                               \
    *obj = (Type){__VA_ARGS__};                                                                                        \
    return obj;                                                                                                        \
  }
#define define_init_new(Type, ...)                                                                                     \
  {                                                                                                                    \
    Type *obj = my_malloc(sizeof(Type));                                                                               \
    Type##_init(obj, ##__VA_ARGS__);                                                                                   \
    return (void *)obj;                                                                                                \
  }
void *my_malloc(size_t size);

#endif // UTILS_H
