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