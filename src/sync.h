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
#include <stdint.h>
#include <stdbool.h>
#include "profile.h"

#ifdef _WIN32
#define WIN32_MEAN_AND_LEAN
#include <windows.h>
typedef CRITICAL_SECTION   os_mutex_t;
typedef CONDITION_VARIABLE os_condvar_t;
#elif defined(__linux__)
#include <pthread.h>
#include <semaphore.h>
typedef pthread_mutex_t os_mutex_t;
typedef pthread_cond_t  os_condvar_t;
#endif

typedef struct {
#ifdef _WIN32
    void *data;
#else
    sem_t data;
#endif
} os_semaphore_t;

typedef struct {
    int count;
    os_mutex_t mutex;
    os_condvar_t cond;
} semaphore_t;

void os_mutex_create(os_mutex_t *mutex);
void os_mutex_delete(os_mutex_t *mutex);
void os_mutex_lock  (os_mutex_t *mutex);
void os_mutex_unlock(os_mutex_t *mutex);

void os_condvar_create(os_condvar_t *condvar);
void os_condvar_delete(os_condvar_t *condvar);
bool os_condvar_wait  (os_condvar_t *condvar, os_mutex_t *mutex, int timeout_ms);
void os_condvar_signal(os_condvar_t *condvar);

void semaphore_create(semaphore_t *sem, int count);
void semaphore_delete(semaphore_t *sem);
bool semaphore_wait  (semaphore_t *sem, int count, int timeout_ms);
void semaphore_signal(semaphore_t *sem, int count);

bool os_semaphore_create(os_semaphore_t *sem, int count, int max);
bool os_semaphore_delete(os_semaphore_t *sem);
bool os_semaphore_wait  (os_semaphore_t *sem);
bool os_semaphore_signal(os_semaphore_t *sem);

profile_results_t sync_profile_results(void);