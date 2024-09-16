#include <stddef.h>

#define COUNTOF(X) (int) (sizeof(X) / sizeof((X)[0]))

char *load_file(const char *file, size_t *size);