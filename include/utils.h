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
    Type##_init(obj, __VA_ARGS__);                                                                                     \
    return obj;                                                                                                        \
  }
void *my_malloc(size_t size);

#define define_list_header(Type)                                                                                       \
  typedef struct Type##List {                                                                                          \
    size_t max_size;                                                                                                   \
    size_t size;                                                                                                       \
    Type *items;                                                                                                       \
  } Type##List;                                                                                                        \
  void Type##List_init(Type##List *list, size_t max_size);                                                             \
  Type##List *Type##List_new(size_t max_size);                                                                         \
  void Type##List_append(Type##List *list, Type item);

#define define_list_source(Type)                                                                                       \
  void Type##List_init(Type##List *list, size_t max_size) {                                                            \
    list->max_size = max_size;                                                                                         \
    list->size = 0;                                                                                                    \
    list->items = my_malloc(sizeof(list->items[0]) * max_size);                                                        \
  }                                                                                                                    \
  Type##List *Type##List_new(size_t max_size) define_init_new(Type##List, max_size);                                   \
  void Type##List_append(Type##List *list, Type item) {                                                                \
    assert((list->size < list->max_size) && "List is full");                                                           \
    list->items[list->size++] = item;                                                                                  \
  }

#endif // UTILS_H
