/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org/>
 */

#include <stdatomic.h>

#ifdef _WIN32
#else
#include <stdlib.h>
#endif

#include "thread.h"

void os_thread_create(os_thread *thread, void *arg, os_threadreturn (*func)(void*))
{
    #if defined(_WIN32)
    os_thread thread_ = CreateThread(NULL, 0, func, arg, 0, NULL);
    if (thread_ == INVALID_HANDLE_VALUE)
        abort();
    *thread = thread_;
    #elif defined(__linux__)
    int ret = pthread_create(thread, NULL, func, arg);
    if (ret) abort();
    #endif
}

os_threadreturn os_thread_join(os_thread thread)
{
    #if defined(_WIN32)
    os_threadreturn result;
    WaitForSingleObject(thread, INFINITE);
    if (!GetExitCodeThread(thread, &result))
        abort();
    CloseHandle(thread);
    return result;
    #elif defined(__linux__)
    os_threadreturn result;
    int ret = pthread_join(thread, &result);
    if (ret) abort();
    return result;
    #else
    (void) thread;
    #endif
}

uint64_t get_thread_id(void)
{
    static _Atomic uint64_t next_id = 1;
    static _Thread_local uint64_t id = 0;
    if (id == 0) id = atomic_fetch_add(&next_id, 1);
    return id;
}