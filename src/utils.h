#include <stddef.h>
#include <stdbool.h>

#define COUNTOF(X) (int) (sizeof(X) / sizeof((X)[0]))

char *load_file(const char *file, size_t *size);

inline bool is_space(char c) { return c == ' ' || c == '\r' || c == '\t' || c == '\n'; }
inline bool is_digit(char c) { return c >= '0' && c <= '9'; }