#include <assert.h>
#include <stdlib.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#define clamp(x, lo, hi) min(max(x, lo), hi)

inline void *my_malloc(size_t size) {
  void *ptr = malloc(size);
  assert((ptr != NULL) && "Failed to allocate memory");
  return ptr;
}
