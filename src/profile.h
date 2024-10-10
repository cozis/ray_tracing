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
#ifndef PROFILE_H
#define PROFILE_H

#include <stdint.h>
#include <x86intrin.h>
#include <stdatomic.h>

typedef struct {
    const char *_Atomic name;
    _Atomic uint64_t elapsed_cycles;
    _Atomic uint64_t call_count;
} profile_t;

typedef struct {
    profile_t *array;
    int count;
} profile_results_t;

#ifdef PROFILE

#define PROFILE_START                                     \
    profile_t *profile__ = &profile_table__[__COUNTER__]; \
    profile__->name = __func__;                            \
    uint64_t profile_start__ = __rdtsc();

#define PROFILE_END \
    atomic_fetch_add(&profile__->elapsed_cycles, __rdtsc() - profile_start__); \
    atomic_fetch_add(&profile__->call_count, 1);

#define PROFILE_GLOBAL_START \
    static profile_t profile_table__[];

#define PROFILE_GLOBAL_END \
    static profile_t profile_table__[__COUNTER__];

#define PROFILE_RESULTS (profile_results_t) {profile_table__, sizeof(profile_table__) / sizeof(profile_table__[0])}

#else
#define PROFILE_START
#define PROFILE_END
#define PROFILE_GLOBAL_START
#define PROFILE_GLOBAL_END
#define PROFILE_RESULTS (profile_results_t) {NULL, 0}
#endif
#endif

void human_readable_time_interval(uint64_t ns, char *dst, size_t max);
void print_profile_results(profile_results_t res_list[], int num_results, long double ns_per_cycle);