#include "utils.h"

void *my_malloc(size_t size) {
  void *ptr = malloc(size);
  assert((ptr != NULL) && "Failed to allocate memory");
  return ptr;
}

void list_init(List *list, size_t max_size) {
  list->max_size = max_size;
  list->size = 0;
  list->items = my_malloc(sizeof(list->items[0]) * max_size);
}
void list_append(List *list, void *obj) {
  assert((list->size <= list->max_size) && ("Max size reached"));
  list->items[list->size++] = obj;
}
void *list_head(List *list) {
  assert((list->size > 0) && ("List is empty"));
  return list->items[list->size - 1];
}
