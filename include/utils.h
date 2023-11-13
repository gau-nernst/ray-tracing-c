#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdlib.h>

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
    Type##_init(obj, __VA_ARGS__);                                                                                     \
    return obj;                                                                                                        \
  }
void *my_malloc(size_t size);

typedef struct List {
  size_t max_size;
  size_t size;
  void **items;
} List;

void list_init(List *list, size_t max_size);
void list_append(List *list, void *obj);
void *list_head(List *list);

#endif // UTILS_H
