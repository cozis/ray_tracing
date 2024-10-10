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
#include <stdlib.h>
#include <assert.h>

#if defined(__linux__)
#include <errno.h>
#endif

#include "sync.h"
#include "clock.h"

//#define SYNC_PRINT_ERRORS
#ifdef SYNC_PRINT_ERRORS
#include <stdio.h>
#include <string.h>
#endif

PROFILE_GLOBAL_START;

void os_mutex_create(os_mutex_t *mutex)
{
    PROFILE_START;

#if defined(_WIN32)
    InitializeCriticalSection(mutex);
#elif defined(__linux__)
    if (pthread_mutex_init(mutex, NULL))
        abort();
#else
    (void) mutex;
#endif

    PROFILE_END;
}

void os_mutex_delete(os_mutex_t *mutex)
{
    PROFILE_START;

#if defined(_WIN32)
    DeleteCriticalSection(mutex);
#elif defined(__linux__)
    if (pthread_mutex_destroy(mutex))
        abort();
#else
    (void) mutex;
#endif

    PROFILE_END;
}

void os_mutex_lock(os_mutex_t *mutex)
{
    PROFILE_START;

#if defined(_WIN32)
    EnterCriticalSection(mutex);
#elif defined(__linux__)
    if (pthread_mutex_lock(mutex))
        abort();
#else
    (void) mutex;
#endif

    PROFILE_END;
}

void os_mutex_unlock(os_mutex_t *mutex)
{
    PROFILE_START;

#if defined(_WIN32)
    LeaveCriticalSection(mutex);
#elif defined(__linux__)
    if (pthread_mutex_unlock(mutex))
        abort();
#else
    (void) mutex;
#endif

    PROFILE_END;
}

void os_condvar_create(os_condvar_t *condvar)
{
    PROFILE_START;

#if defined(_WIN32)
    InitializeConditionVariable(condvar);
#elif defined(__linux__)
    if (pthread_cond_init(condvar, NULL))
        abort();
#else
    (void) condvar;
#endif

    PROFILE_END;
}

void os_condvar_delete(os_condvar_t *condvar)
{
    PROFILE_START;

#if defined(__linux__)
    if (pthread_cond_destroy(condvar))
        abort();
#else
    (void) condvar;
#endif

    PROFILE_END;
}

bool os_condvar_wait(os_condvar_t *condvar, os_mutex_t *mutex, int timeout_ms)
{
    PROFILE_START;

#if defined(_WIN32)

    DWORD timeout = INFINITE;
    if (timeout_ms >= 0) timeout = timeout_ms;
    if (!SleepConditionVariableCS(condvar, mutex, timeout)) {
        if (GetLastError() == ERROR_TIMEOUT) {
            PROFILE_END;
            return false;
        }
        abort();
    }

    PROFILE_END;
    return true;

#elif defined(__linux__)

    int err;
    if (timeout_ms < 0)
        err = pthread_cond_wait(condvar, mutex);
    else {
        uint64_t wakeup_ms = (uint64_t) timeout_ms + get_absolute_time_us() / 1000;
        struct timespec abstime = {
            .tv_sec = wakeup_ms / 1000,
            .tv_nsec = (wakeup_ms % 1000) * 1000000,
        };
        err = pthread_cond_timedwait(condvar, mutex, &abstime);
    }
    if (err) {
        if (err == ETIMEDOUT) {
            PROFILE_END;
            return false;
        }
#ifdef SYNC_PRINT_ERRORS
        fprintf(stderr, "ERROR!! pthread_cond_wait/timedwait: %s\n", strerror(err));
#endif
        abort();
    }

    PROFILE_END;
    return true;

#else
    PROFILE_END;
    (void) condvar;
#endif
}

void os_condvar_signal(os_condvar_t *condvar)
{
    PROFILE_START;

#if defined(_WIN32)
    WakeConditionVariable(condvar);
#elif defined(__linux__)
    if (pthread_cond_signal(condvar))
        abort();
#else
    (void) condvar;
#endif

    PROFILE_END;
}

void semaphore_create(semaphore_t *sem, int count)
{
    PROFILE_START;

    sem->count = count;
    os_mutex_create(&sem->mutex);
    os_condvar_create(&sem->cond);

    PROFILE_END;
}

void semaphore_delete(semaphore_t *sem)
{
    PROFILE_START;

    os_mutex_delete(&sem->mutex);
    os_condvar_delete(&sem->cond);

    PROFILE_END;
}

bool semaphore_wait(semaphore_t *sem, int count, int timeout_ms)
{
    PROFILE_START;

    assert(count > 0);

    uint64_t start_time_ms = get_relative_time_ns() / 1000000;

    os_mutex_lock(&sem->mutex);
    while (sem->count < count) {

        uint64_t current_time_ms = get_relative_time_ns() / 1000000;
        
        int remaining_ms = -1;
        if (timeout_ms >= 0)
            remaining_ms = timeout_ms - (int) (current_time_ms - start_time_ms);
        
        if (!os_condvar_wait(&sem->cond, &sem->mutex, remaining_ms)) {
            os_mutex_unlock(&sem->mutex);
            PROFILE_END;
            return false;
        }
    }
    sem->count -= count;
    os_mutex_unlock(&sem->mutex);

    PROFILE_END;
    return true;
}

void semaphore_signal(semaphore_t *sem, int count)
{
    PROFILE_START;

    assert(count > 0);

    os_mutex_lock(&sem->mutex);
    sem->count += count;
    if (sem->count > 0)
        os_condvar_signal(&sem->cond);
    os_mutex_unlock(&sem->mutex);

    PROFILE_END;
}

bool os_semaphore_create(os_semaphore_t *sem, int count, int max)
{
    PROFILE_START;

    int ok;

    #ifdef _WIN32
    SECURITY_ATTRIBUTES *attr = NULL; // Default
    const char *name = NULL; // No name
    void *handle = CreateSemaphoreA(attr, count, max, name);
    if (handle == NULL) {
        PROFILE_END;
        return false;
    }
    sem->data = handle;
    ok = 1;
    #else
    (void) max; // POSIX doesn't use this
    ok = sem_init(&sem->data, 0, count) == 0;
    #endif

    PROFILE_END;
    return ok;
}

bool os_semaphore_delete(os_semaphore_t *sem)
{
    PROFILE_START;

    int ok;

    #ifdef _WIN32
    CloseHandle(sem->data);
    ok = 1;
    #else
    ok = sem_destroy(&sem->data) == 0;
    #endif

    PROFILE_END;
    return ok;
}

bool os_semaphore_wait(os_semaphore_t *sem)
{
    PROFILE_START;

    int ok;

    #ifdef _WIN32
    ok = WaitForSingleObject(sem->data, INFINITE) == WAIT_OBJECT_0;
    #else
    ok = sem_wait(&sem->data) == 0;
    #endif

    PROFILE_END;
    return ok;
}

bool os_semaphore_signal(os_semaphore_t *sem)
{
    PROFILE_START;

    int ok;

    #ifdef _WIN32
    ok = ReleaseSemaphore(sem->data, 1, NULL);
    #else
    ok = sem_post(&sem->data) == 0;
    #endif

    PROFILE_END;
    return ok;
}

PROFILE_GLOBAL_END;

profile_results_t sync_profile_results(void)
{
    return PROFILE_RESULTS;
}