#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

char *load_file(const char *file, size_t *size)
{
    FILE *stream = fopen(file, "rb");
    if (stream == NULL) return NULL;

    fseek(stream, 0, SEEK_END);
    long size2 = ftell(stream);
    fseek(stream, 0, SEEK_SET);

    char *dst = (char*) malloc(size2+1);
    if (dst == NULL) {
        fclose(stream);
        return NULL;
    }

    fread(dst, 1, size2, stream);
    if (ferror(stream)) {
        free(dst);
        fclose(stream);
        return NULL;
    }
    dst[size2] = '\0';

    fclose(stream);
    if (size) *size = size2;
    return dst;
}
