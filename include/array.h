#include "my_utils.h"
#include <assert.h>
#include <stdlib.h>

#define define_array_typedef(Type)                                                                                     \
  typedef struct Type##Array {                                                                                         \
    size_t max_size;                                                                                                   \
    size_t size;                                                                                                       \
    Type *items;                                                                                                       \
  } Type##Array;

#define define_array_new(Type)                                                                                         \
  Type##Array Type##Array_new(size_t max_size) {                                                                       \
    return (Type##Array){max_size, 0, my_malloc(sizeof(Type) * max_size)};                                             \
  }

#define define_array_append(Type)                                                                                      \
  Type *Type##Array_append(Type##Array *array, Type object) {                                                          \
    assert((array->size < array->max_size) && "Max size reached");                                                     \
    array->items[array->size++] = object;                                                                              \
    return array->items + array->size - 1;                                                                             \
  }

#define define_array_next(Type)                                                                                        \
  Type *Type##Array_next(Type##Array *array) {                                                                         \
    assert((array->size < array->max_size) && "Max size reached");                                                     \
    array->size++;                                                                                                     \
    return array->items + array->size - 1;                                                                             \
  }

#define define_array_source(Type)                                                                                      \
  define_array_new(Type);                                                                                              \
  define_array_append(Type);                                                                                           \
  define_array_next(Type);

#define define_array_header(Type)                                                                                      \
  define_array_typedef(Type);                                                                                          \
  Type##Array Type##Array_new(size_t max_size);                                                                        \
  Type *Type##Array_append(Type##Array *array, Type object);                                                           \
  Type *Type##Array_next(Type##Array *array)
