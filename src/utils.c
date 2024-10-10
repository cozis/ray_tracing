/*
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
*/
#include <stdio.h>
#include <stdint.h>
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

static _Thread_local uint64_t wyhash64_x = 0;

static uint64_t wyhash64(void) {
	wyhash64_x += 0x60bee2bee120fc15;
	__uint128_t tmp;
	tmp = (__uint128_t) wyhash64_x * 0xa3b195354a39b70d;
	uint64_t m1 = (tmp >> 64) ^ tmp;
	tmp = (__uint128_t)m1 * 0x1b03738712fad5c9;
	uint64_t m2 = (tmp >> 64) ^ tmp;
	return m2;
}

float random_float(void)
{
	return (float) wyhash64() / UINT64_MAX;
}
