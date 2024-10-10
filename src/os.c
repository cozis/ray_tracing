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
#define _GNU_SOURCE
#define _POSIX_C_SOURCE 1999309L

#include <stdlib.h>
#include <assert.h>
#include <stdatomic.h>

#ifdef _WIN32
#define WIN32_MEAN_AND_LEAN
#include <windows.h>
#endif

#ifdef __linux__
#include <time.h>
#include <errno.h>
#include <unistd.h>
#endif

//#define SYNC_PRINT_ERRORS
#ifdef SYNC_PRINT_ERRORS
#include <stdio.h>
#include <string.h>
#endif

#include "os.h"

/*
 * Returns the current absolute time in microsecods
 * TODO: Specify since when the time is calculated
 */
uint64_t get_absolute_time_us(void)
{
#ifdef _WIN32
	FILETIME filetime;
	GetSystemTimePreciseAsFileTime(&filetime);
	uint64_t time = (uint64_t) filetime.dwLowDateTime | ((uint64_t) filetime.dwHighDateTime << 32);
	time /= 10;
	return time;
#else
	struct timespec buffer;
	if (clock_gettime(CLOCK_REALTIME, &buffer))
		abort();
	uint64_t time = buffer.tv_sec * 1000000 + buffer.tv_nsec / 1000;
	return time;
#endif
}

/*
 * Returns the current time in nanoseconds since
 * an unspecified point in time.
 */
uint64_t get_relative_time_ns(void)
{
#ifdef _WIN32
	{
		int64_t count;
		int64_t freq;
		int ok;
		
		ok = QueryPerformanceCounter((LARGE_INTEGER*) &count);
		if (!ok) abort();

		ok = QueryPerformanceFrequency((LARGE_INTEGER*) &freq);
		if (!ok) abort();

		uint64_t res = 1000000000 * (double) count / freq;
		return res;
	}
#else
	{
		struct timespec time;

		if (clock_gettime(CLOCK_REALTIME, &time))
			abort();

		uint64_t res;

		uint64_t sec = time.tv_sec;
		if (sec > UINT64_MAX / 1000000000)
			abort();
		res = sec * 1000000000;
		
		uint64_t nsec = time.tv_nsec;
		if (res > UINT64_MAX - nsec)
			abort();
		res += nsec;

		return res;
	}
#endif
}

void sleep_ms(float ms)
{
#ifdef _WIN32
	Sleep(ms);
#else
	usleep(ms*1000);
#endif
}

void os_mutex_create(os_mutex_t *mutex)
{

#if defined(_WIN32)
	InitializeCriticalSection(mutex);
#elif defined(__linux__)
	if (pthread_mutex_init(mutex, NULL))
		abort();
#else
	(void) mutex;
#endif

}

void os_mutex_delete(os_mutex_t *mutex)
{

#if defined(_WIN32)
	DeleteCriticalSection(mutex);
#elif defined(__linux__)
	if (pthread_mutex_destroy(mutex))
		abort();
#else
	(void) mutex;
#endif

}

void os_mutex_lock(os_mutex_t *mutex)
{

#if defined(_WIN32)
	EnterCriticalSection(mutex);
#elif defined(__linux__)
	if (pthread_mutex_lock(mutex))
		abort();
#else
	(void) mutex;
#endif

}

void os_mutex_unlock(os_mutex_t *mutex)
{

#if defined(_WIN32)
	LeaveCriticalSection(mutex);
#elif defined(__linux__)
	if (pthread_mutex_unlock(mutex))
		abort();
#else
	(void) mutex;
#endif

}

void os_condvar_create(os_condvar_t *condvar)
{

#if defined(_WIN32)
	InitializeConditionVariable(condvar);
#elif defined(__linux__)
	if (pthread_cond_init(condvar, NULL))
		abort();
#else
	(void) condvar;
#endif

}

void os_condvar_delete(os_condvar_t *condvar)
{

#if defined(__linux__)
	if (pthread_cond_destroy(condvar))
		abort();
#else
	(void) condvar;
#endif

}

bool os_condvar_wait(os_condvar_t *condvar, os_mutex_t *mutex, int timeout_ms)
{

#if defined(_WIN32)

	DWORD timeout = INFINITE;
	if (timeout_ms >= 0) timeout = timeout_ms;
	if (!SleepConditionVariableCS(condvar, mutex, timeout)) {
		if (GetLastError() == ERROR_TIMEOUT) {
			return false;
		}
		abort();
	}

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
			return false;
		}
	#ifdef SYNC_PRINT_ERRORS
		fprintf(stderr, "ERROR!! pthread_cond_wait/timedwait: %s\n", strerror(err));
	#endif
		abort();
	}

	return true;

#else
	(void) condvar;
#endif
}

void os_condvar_signal(os_condvar_t *condvar)
{

#if defined(_WIN32)
	WakeConditionVariable(condvar);
#elif defined(__linux__)
	if (pthread_cond_signal(condvar))
		abort();
#else
    (void) condvar;
#endif

}

void semaphore_create(semaphore_t *sem, int count)
{
	sem->count = count;
	os_mutex_create(&sem->mutex);
	os_condvar_create(&sem->cond);

}

void semaphore_delete(semaphore_t *sem)
{
	os_mutex_delete(&sem->mutex);
	os_condvar_delete(&sem->cond);
}

bool semaphore_wait(semaphore_t *sem, int count, int timeout_ms)
{
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
			return false;
		}
	}
	sem->count -= count;
	os_mutex_unlock(&sem->mutex);

	return true;
}

void semaphore_signal(semaphore_t *sem, int count)
{
	assert(count > 0);

	os_mutex_lock(&sem->mutex);
	sem->count += count;
	if (sem->count > 0)
		os_condvar_signal(&sem->cond);
	os_mutex_unlock(&sem->mutex);
}

bool os_semaphore_create(os_semaphore_t *sem, int count, int max)
{
	int ok;

#ifdef _WIN32
	SECURITY_ATTRIBUTES *attr = NULL; // Default
	const char *name = NULL; // No name
	void *handle = CreateSemaphoreA(attr, count, max, name);
	if (handle == NULL)
		return false;
	sem->data = handle;
	ok = 1;
#else
	(void) max; // POSIX doesn't use this
	ok = sem_init(&sem->data, 0, count) == 0;
#endif

	return ok;
}

bool os_semaphore_delete(os_semaphore_t *sem)
{
	int ok;

#ifdef _WIN32
	CloseHandle(sem->data);
	ok = 1;
#else
	ok = sem_destroy(&sem->data) == 0;
#endif

	return ok;
}

bool os_semaphore_wait(os_semaphore_t *sem)
{
	int ok;

#ifdef _WIN32
	ok = WaitForSingleObject(sem->data, INFINITE) == WAIT_OBJECT_0;
#else
	ok = sem_wait(&sem->data) == 0;
#endif

	return ok;
}

bool os_semaphore_signal(os_semaphore_t *sem)
{
	int ok;

#ifdef _WIN32
	ok = ReleaseSemaphore(sem->data, 1, NULL);
#else
	ok = sem_post(&sem->data) == 0;
#endif

	return ok;
}

void os_thread_create(os_thread *thread, void *arg, os_threadreturn (*func)(void*))
{
#ifdef _WIN32
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
#ifdef _WIN32
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