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
#include "profile.h"

void human_readable_time_interval(uint64_t ns, char *dst, size_t max)
{
    if (ns > 1000000000)
        snprintf(dst, max, "%.1Lf s", (long double) ns / 1000000000);
    else if (ns > 1000000)
        snprintf(dst, max, "%.1Lf ms", (long double) ns / 1000000);
    else if (ns > 1000)
        snprintf(dst, max, "%.1Lf us", (long double) ns / 1000);
    else
        snprintf(dst, max, "%.1Lf ns", (long double) ns);
}

static int
sort_callback(const void *a, const void *b)
{
    const profile_t *p1 = a;
    const profile_t *p2 = b;
    uint64_t n1 = atomic_load(&p1->elapsed_cycles);
    uint64_t n2 = atomic_load(&p2->elapsed_cycles);
    if (n1 > n2) return +1;
    if (n1 < n2) return -1;
    return 0;
}

void print_profile_results(profile_results_t res_list[], int num_results, long double ns_per_cycle)
{
    int entry_count = 0;
    for (int i = 0; i < num_results; i++)
        entry_count += res_list[i].count;

    profile_t  entries_[128];
    profile_t *entries = entries_;

    if (entry_count >= (int) (sizeof(entries_)/sizeof(profile_t))) {
        entries = malloc(entry_count * sizeof(profile_t));
        if (entries == NULL) abort();
    }

    for (int i = 0, k = 0; i < num_results; i++)
        for (int j = 0; j < res_list[i].count; j++)
            entries[k++] = res_list[i].array[j];

    qsort(entries, entry_count, sizeof(profile_t), sort_callback);

    if (entry_count > 0)
        fprintf(stderr, "| %-30s | %-10s | %-10s | %-10s |\n",
            "Name", "Calls", "Total", "Latency");


    for (int i = 0; i < entry_count; i++) {
        if (!entries[i].name) continue;
        int call_count = entries[i].call_count;
        uint64_t elapsed_cycles = entries[i].elapsed_cycles;

        long double total_ns = ns_per_cycle * elapsed_cycles;
        long double latency_ns = total_ns / call_count;

        char total[128];
        human_readable_time_interval(total_ns, total, sizeof(total));

        char latency[128];
        human_readable_time_interval(latency_ns, latency, sizeof(latency));

        fprintf(stderr, "| %-30s | %-10d | %-10s | %-10s |\n",
            entries[i].name, call_count, total, latency);
    }
}