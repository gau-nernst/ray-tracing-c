#include "utils.h"

void *my_malloc(size_t size) {
  void *ptr = malloc(size);
  assert((ptr != NULL) && "Failed to allocate memory");
  return ptr;
}
