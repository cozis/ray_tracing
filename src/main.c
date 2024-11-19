/*
Copyright 2024 Francesco Cozzuto

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the “Software”), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdatomic.h>
#include <x86intrin.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "os.h"
#include "utils.h"
#include "camera.h"
#include "scene.h"
#include "gpu_and_windowing.h"

typedef struct {
	int column_i;
} WorkerConfig;

typedef struct {
	float   depth;
	Vector3 color;
	Vector3 normal;
	float   roughness;
} PixelInfo;

typedef void (*FilterFunc)(Vector3*, Vector3*, float*, float*, Vector3*, int, int, int, int, int);
/////////////////////////////////////////////////////////////////////////////
/// GLOBAL VARIABLES                                                      ///
/////////////////////////////////////////////////////////////////////////////

#define MAX_COLUMNS 32
#define INIT_SMOOTH_LIMIT 10
#define SMOOTH_RATIO 1.5

// Parameters. These are set at startup and are
// considered constant after that.
int num_columns;

// The scene and background being rendered.
Scene   scene;
Cubemap skybox;

int color_buffer_frames = 0;
Vector3 *color_buffer;
float   *depth_buffer;
Vector3 *normal_buffer;
Vector3 *filter_buffer_A;
Vector3 *filter_buffer_B;
float   *roughness_buffer;

// Size of the accumulation and frame buffers
int frame_w = 0;
int frame_h = 0;

// This guards the critical section around the accumulation buffer.
os_mutex_t frame_mutex;

int completed = 0;
os_condvar_t start_work;
os_condvar_t completed_work;

uint64_t loop_cycles = 0;
uint64_t loop_count = 0;

uint64_t global_frame_index = 0;

uint64_t global_cycle_pixel_sum = 0;
uint64_t global_pixel_count = 0;

uint64_t smooth_index = 0;
uint64_t smooth_limit = INIT_SMOOTH_LIMIT;

/////////////////////////////////////////////////////////////////////////////
/// FUNCTION PROTOTYPES                                                   ///
/////////////////////////////////////////////////////////////////////////////

void start_filtering_workers(void);
void stop_filtering_workers(void);
void filter(Vector3 *input, Vector3 *output, FilterFunc func);

void    start_ray_tracing_workers(void);
void    stop_ray_tracing_workers(void);

bool    quitting(void);
void    screenshot(void);
void    parse_arguments_or_exit(int argc, char **argv, int *num_columns, char **scene_file);

PixelInfo pixel(float x, float y, float aspect_ratio);
void      update_frame(void);
void      invalidate_accumulation(void);

os_threadreturn worker(void *arg);

/////////////////////////////////////////////////////////////////////////////
/// IMPLEMENTATION                                                        ///
/////////////////////////////////////////////////////////////////////////////

void invalidate_accumulation(void)
{
	os_mutex_lock(&frame_mutex);
	memset(color_buffer,  0, sizeof(Vector3) * frame_w * frame_h);
	memset(depth_buffer,  0, sizeof(float)   * frame_w * frame_h);
	memset(normal_buffer, 0, sizeof(Vector3) * frame_w * frame_h);
	memset(roughness_buffer, 0, sizeof(float) * frame_w * frame_h);
	color_buffer_frames = 0;
	smooth_index = 0;
	smooth_limit = INIT_SMOOTH_LIMIT;
	os_mutex_unlock(&frame_mutex);
}

Vector3 fresnel_schlick(float u, Vector3 f0)
{
	return combine(f0, combine(vec_from_scalar(1.0), f0, 1, -1), 1, pow(1.0 - u, 5.0));
}

Vector3 cosine_sample(Vector3 n)
{
	float eps = 0.001;
	Vector3 u = {1, 0, 0};
	if (absf(n.y) < 1 - eps)
		u = normalize(cross(n, (Vector3) {0, 1, 0}));
	Vector3 v = cross(n, u);
	float phi = 2.0 * M_PI * random_float();
	float ay = sqrt(random_float());
	float ax = sqrt(1.0 - ay*ay);

	Vector3 result;
	result.x = ax * (cos(phi) * u.x + sin(phi) * v.x) + ay * n.x;
	result.y = ax * (cos(phi) * u.y + sin(phi) * v.y) + ay * n.y;
	result.z = ax * (cos(phi) * u.z + sin(phi) * v.z) + ay * n.z;
	return result;
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

uint64_t rand64(void)
{
	return wyhash64();
}

PixelInfo pixel_inner(Ray in_ray)
{
	// Find a light source. This is kind of lazy as we should
	// sample every light source in the scene.
	float light_sample_weight = 0.05;
	Vector3 weighted_light_emission = {0, 0, 0};
	Vector3 light_emission = {0, 0, 0};
	Object *light_object = NULL;
	Vector3 light_origin;
	for (int i = 0; i < scene.num_objects; i++) {
		Material material = scene.objects[i].material;
		if (material.emission_power > 0) {
			weighted_light_emission.x += material.emission_color.x * material.emission_power * light_sample_weight;
			weighted_light_emission.y += material.emission_color.y * material.emission_power * light_sample_weight;
			weighted_light_emission.z += material.emission_color.z * material.emission_power * light_sample_weight;
			light_emission.x += material.emission_color.x * material.emission_power;
			light_emission.y += material.emission_color.y * material.emission_power;
			light_emission.z += material.emission_color.z * material.emission_power;
			light_object = &scene.objects[i];
			light_origin = origin_of(scene.objects[i]);
			break;
		}
	}

	// Keep track of how much of the light ray has been
	// absorbed while bouncing around
	Vector3 contrib = {1, 1, 1};

	// Keep track of the final luminosity
	Vector3 result = {0, 0, 0};

	// Maximum number of bounces of the ray
	int bounces = 6;

	PixelInfo info;
	info.depth = 1000000;
	info.color = (Vector3) {0, 0, 0};
	info.normal = (Vector3) {-in_ray.direction.x, -in_ray.direction.y, -in_ray.direction.z};
	info.roughness = 0;

	for (int i = 0; i < bounces; i++) {

		// Find the next collision
		HitInfo hit = trace_ray(in_ray, &scene);
		if (hit.object == -1) {
			// The ray flew straight out of the scene!
			//
			// The sky is sampled here. You can change the sky color here
			// if you want:
			//     Vector3 sky_color = {0.6, 0.7, 0.9};
			     Vector3 sky_color = {0, 0, 0};
			//     Vector3 sky_color = {1, 1, 1};
			//Vector3 sky_color = sample_cubemap(&skybox, normalize(in_ray.direction));
			result = combine(result, mulv(sky_color, contrib), 1, 1);
			break;
		}

		Material material = scene.objects[hit.object].material;

		if (i == 0) {
			info.depth = hit.distance;
			info.normal = hit.normal;
			info.roughness = material.roughness;
		}

		uint64_t rand_bucket_0 = rand64();

		Vector3 v = in_ray.direction;
		Vector3 n = hit.normal;
		Vector3 o = hit.point;

		{
			float norm2 = v.x*v.x + v.y*v.y + v.z*v.z;
			float norm = sqrt(norm2);
			v.x /= norm;
			v.y /= norm;
			v.z /= norm;
		}

		uint64_t rand_0 = (rand_bucket_0 >>  0) & 0xFFFFF;
		uint64_t rand_1 = (rand_bucket_0 >> 20) & 0xFFFFF;
		uint64_t rand_2 = (rand_bucket_0 >> 40) & 0xFFFFF;

		uint64_t rand_bucket_1 = rand64();

		result.x += contrib.x * material.emission_color.x * material.emission_power;
		result.y += contrib.y * material.emission_color.y * material.emission_power;
		result.z += contrib.z * material.emission_color.z * material.emission_power;

		float NoV = -(n.x * v.x + n.y * v.y + n.z * v.z);
		if (NoV < 0)
			NoV = 0;
		else {
			if (NoV > 1)
				NoV = 1;
		}

		float f0 = 0.16 * material.reflectance * material.reflectance * (1 - material.metallic) + avgv(material.albedo) * material.metallic;
		float F = f0 + (1 - f0) * pow(1 - NoV, 5);

		Vector3 rand_dir;
		float  rand_dir_factor;
		{
#if 0
			rand_dir_factor = (float) 2 / 0xFFFFF;

			rand_dir.x = (float) rand_0 - 0.5 * 0xFFFFF;
			rand_dir.y = (float) rand_1 - 0.5 * 0xFFFFF;
			rand_dir.z = (float) rand_2 - 0.5 * 0xFFFFF;

			if (rand_dir.x * n.x + rand_dir.y * n.y + rand_dir.z * n.z < 0)
				rand_dir_factor = -rand_dir_factor;
#else
			Vector3 dir_to_light;
			dir_to_light.x = light_origin.x - o.x;
			dir_to_light.y = light_origin.y - o.y;
			dir_to_light.z = light_origin.z - o.z;

			rand_dir = cosine_sample(dir_to_light);
			float p = dotv(rand_dir, n);
			assert(p > 0);
			rand_dir_factor = 1;

			contrib.x /= p;
			contrib.y /= p;
			contrib.z /= p;

#endif
		}

		Vector3 out_dir;
		uint64_t rand_6 = rand64();
		if (material.metallic > 0.001 || rand_6 <= F * UINT64_MAX) {

			// Specular ray

			Vector3 reflect_dir;
			float tmp = 2 * NoV;
			reflect_dir.x = tmp * n.x + v.x;
			reflect_dir.y = tmp * n.y + v.y;
			reflect_dir.z = tmp * n.z + v.z;

			rand_dir_factor *= material.roughness;
			out_dir.x = rand_dir_factor * rand_dir.x + reflect_dir.x;
			out_dir.y = rand_dir_factor * rand_dir.y + reflect_dir.y;
			out_dir.z = rand_dir_factor * rand_dir.z + reflect_dir.z;

		} else {

			// Diffuse ray

			if (false) {

				Vector3 dir_to_light;
				dir_to_light.x = light_origin.x - o.x;
				dir_to_light.y = light_origin.y - o.y;
				dir_to_light.z = light_origin.z - o.z;

				Vector3 sample_dir = cosine_sample(dir_to_light);

				float p = dotv(sample_dir, dir_to_light);
				assert(p < 1.1);

				out_dir.x = sample_dir.x;
				out_dir.y = sample_dir.y;
				out_dir.z = sample_dir.z;

				contrib.x = (1 - material.metallic) * contrib.x * material.albedo.x / p;
				contrib.y = (1 - material.metallic) * contrib.y * material.albedo.y / p;
				contrib.z = (1 - material.metallic) * contrib.z * material.albedo.z / p;

			} else {

				out_dir.x = rand_dir_factor * rand_dir.x;
				out_dir.y = rand_dir_factor * rand_dir.y;
				out_dir.z = rand_dir_factor * rand_dir.z;

				contrib.x = (1 - material.metallic) * contrib.x * material.albedo.x;
				contrib.y = (1 - material.metallic) * contrib.y * material.albedo.y;
				contrib.z = (1 - material.metallic) * contrib.z * material.albedo.z;
			}
		}

		Ray out_ray;
		out_ray.direction = out_dir;
		out_ray.origin.x = hit.point.x + out_dir.x * 0.001;
		out_ray.origin.y = hit.point.y + out_dir.y * 0.001;
		out_ray.origin.z = hit.point.z + out_dir.z * 0.001;

		if (false) {

			Vector3 dir_to_light;
			dir_to_light.x = light_origin.x - o.x;
			dir_to_light.y = light_origin.y - o.y;
			dir_to_light.z = light_origin.z - o.z;

			if (dir_to_light.x*n.x +
				dir_to_light.y*n.y +
				dir_to_light.z*n.z > 0) {
#if 1
				Vector3 sample_dir = cosine_sample(dir_to_light);

				float p = dotv(sample_dir, dir_to_light);
				assert(p < 1.1);

				Ray sample_ray;
				float eps = 0.001;
				sample_ray.direction = sample_dir;
				sample_ray.origin.x = o.x + eps * sample_dir.x;
				sample_ray.origin.y = o.y + eps * sample_dir.y;
				sample_ray.origin.z = o.z + eps * sample_dir.z;

				HitInfo hit2 = trace_ray(sample_ray, &scene);
				if (hit2.object == light_object - scene.objects) {

					result.x += contrib.x * light_emission.x / p;
					result.y += contrib.y * light_emission.y / p;
					result.z += contrib.z * light_emission.z / p;
/*
					contrib.x *= 1 - p;
					contrib.y *= 1 - p;
					contrib.z *= 1 - p;
*/
				}
#else
				float spread = 0;

				uint64_t rand_3 = (rand_bucket_1 >>  0) & 0xFFFFF;
				uint64_t rand_4 = (rand_bucket_1 >> 20) & 0xFFFFF;
				uint64_t rand_5 = (rand_bucket_1 >> 40) & 0xFFFFF;

				Vector3 rand_dir;
				rand_dir.x = (float) rand_3 - 0.5 / 0xFFFFF;
				rand_dir.y = (float) rand_4 - 0.5 / 0xFFFFF;
				rand_dir.z = (float) rand_5 - 0.5 / 0xFFFFF;

				float rand_dir_factor = spread * 2 / 0xFFFFF;
				float sign = 1;

				float dot = rand_dir.x * n.x + rand_dir.y * n.y + rand_dir.z * n.z;
				if (dot < 0) sign = -1;

				Vector3 sample_dir;
				sample_dir.x = sign * rand_dir_factor * rand_dir.x + dir_to_light.x;
				sample_dir.y = sign * rand_dir_factor * rand_dir.y + dir_to_light.y;
				sample_dir.z = sign * rand_dir_factor * rand_dir.z + dir_to_light.z;

				Ray sample_ray;
				float eps = 0.001;
				sample_ray.direction = sample_dir;
				sample_ray.origin.x = o.x + eps * sample_dir.x;
				sample_ray.origin.y = o.y + eps * sample_dir.y;
				sample_ray.origin.z = o.z + eps * sample_dir.z;

				HitInfo hit2 = trace_ray(sample_ray, &scene);
				if (hit2.object == light_object - scene.objects) {
					
					result.x += contrib.x * weighted_light_emission.x / light_distance2;
					result.y += contrib.y * weighted_light_emission.y / light_distance2;
					result.z += contrib.z * weighted_light_emission.z / light_distance2;

					contrib.x *= 1 - light_sample_weight;
					contrib.y *= 1 - light_sample_weight;
					contrib.z *= 1 - light_sample_weight;
				}
#endif
			}
		}

		in_ray = out_ray;
	}

	info.color = result;
	return info;
}

PixelInfo pixel(float x, float y, float aspect_ratio)
{
	assert(!isnan(aspect_ratio));

	Ray in_ray = ray_through_screen_at(x, y, aspect_ratio);
	assert(!isnanv(in_ray.direction));

	return pixel_inner(in_ray);
}

os_threadreturn ray_tracing_worker(void *arg)
{
	// The actual pixels
	Vector3 *local_color_buffer = NULL;
	float   *local_depth_buffer = NULL;
	Vector3 *local_normal_buffer = NULL;
	float   *local_roughness_buffer = NULL;

	// The screen is divided in "num_columns" columns
	int column_i = (int) arg;
	int column_w;

	// Workers need to know the frame size while evaluating pixel
	// values. Since the frame size may change at any time, threads
	// cache their value.
	int cached_frame_w;
	int cached_frame_h;

	uint64_t local_frame_index = 0;

	os_mutex_lock(&frame_mutex);
	for (;;) {

		while (local_frame_index == global_frame_index && !quitting())
			os_condvar_wait(&start_work, &frame_mutex, -1);
		local_frame_index = global_frame_index;

		if (quitting()) {
			completed++;
			os_condvar_signal(&completed_work);
			break;
		}

		// Cache data and check if we need to resize the column buffer
		bool resize = false;
		if (local_color_buffer == NULL)
			resize = true;
		column_w = frame_w / num_columns;
		cached_frame_w = frame_w;
		cached_frame_h = frame_h;

		os_mutex_unlock(&frame_mutex);

		// We need to resize
		if (resize) {

			free(local_color_buffer);
			local_color_buffer = malloc(sizeof(Vector3) * column_w * cached_frame_h);
			if (!local_color_buffer) abort();

			free(local_depth_buffer);
			local_depth_buffer = malloc(sizeof(float) * column_w * cached_frame_h);
			if (!local_depth_buffer) abort();

			free(local_normal_buffer);
			local_normal_buffer = malloc(sizeof(Vector3) * column_w * cached_frame_h);
			if (!local_normal_buffer) abort();

			free(local_roughness_buffer);
			local_roughness_buffer = malloc(sizeof(float) * column_w * cached_frame_h);
			if (!local_roughness_buffer) abort();
		}

		int column_x = column_w * column_i;
		float aspect_ratio = (float) frame_w / frame_h;

		uint64_t local_cycle_pixel_sum = 0;
		uint64_t local_pixel_count = 0;

		// Iterate over each low resolution pixel
		for (int j = 0; j < frame_h; j++) {
			for (int i = 0; i < column_w; i++) {

				float u = (float) (column_x + i) / (frame_w - 1);
				float v = (float) j  / (frame_h - 1);
				u = 1 - u;
				v = 1 - v;

				uint64_t start_time = __rdtsc();
				PixelInfo info = pixel(u, v, aspect_ratio);
				local_cycle_pixel_sum += __rdtsc() - start_time;
				local_pixel_count++;

				local_color_buffer[j * column_w + i] = info.color;
				local_depth_buffer[j * column_w + i] = info.depth;
				local_normal_buffer[j * column_w + i] = info.normal;
				local_roughness_buffer[j * column_w + i] = info.roughness;
			}
		}

		// Now we try publishing the changes
		os_mutex_lock(&frame_mutex);

		global_cycle_pixel_sum += local_cycle_pixel_sum;
		global_pixel_count += local_pixel_count;

		for (int j = 0; j < frame_h; j++)
			for (int i = 0; i < column_w; i++) {

				int src_index = j * column_w + i;
				int dst_index = j * frame_w + (i + column_x);

				color_buffer[dst_index].x = color_buffer[dst_index].x * color_buffer_frames / (color_buffer_frames + 1) + local_color_buffer[src_index].x / (color_buffer_frames + 1);
				color_buffer[dst_index].y = color_buffer[dst_index].y * color_buffer_frames / (color_buffer_frames + 1) + local_color_buffer[src_index].y / (color_buffer_frames + 1);
				color_buffer[dst_index].z = color_buffer[dst_index].z * color_buffer_frames / (color_buffer_frames + 1) + local_color_buffer[src_index].z / (color_buffer_frames + 1);

				depth_buffer[dst_index] = local_depth_buffer[src_index];
				normal_buffer[dst_index] = local_normal_buffer[src_index];
				roughness_buffer[dst_index] = local_roughness_buffer[src_index];
			}

		completed++;
		os_condvar_signal(&completed_work);
	}
	os_mutex_unlock(&frame_mutex);
}

void realloc_frame_buffer(void)
{
	frame_w = get_screen_w();
	frame_h = get_screen_h();

	free(color_buffer);
	color_buffer = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!color_buffer) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	free(depth_buffer);
	depth_buffer = malloc(sizeof(float) * frame_w * frame_h);
	if (!depth_buffer) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	free(normal_buffer);
	normal_buffer = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!normal_buffer) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	free(roughness_buffer);
	roughness_buffer = malloc(sizeof(float) * frame_w * frame_h);
	if (!roughness_buffer) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	free(filter_buffer_A);
	filter_buffer_A = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!filter_buffer_A) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	free(filter_buffer_B);
	filter_buffer_B = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!filter_buffer_B) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	memset(color_buffer,  0, sizeof(Vector3) * frame_w * frame_h);
	memset(depth_buffer,  0, sizeof(float)   * frame_w * frame_h);
	memset(normal_buffer, 0, sizeof(Vector3) * frame_w * frame_h);
	memset(roughness_buffer, 0, sizeof(float) * frame_w * frame_h);

}

bool frame_buffer_size_doesnt_match_window(void)
{
	return frame_w != get_screen_w() || frame_h != get_screen_h();
}

static int compare_vector_luminosity(const void *a, const void *b)
{
	Vector3 u = *(Vector3*) a;
	Vector3 v = *(Vector3*) b;
	float x = avgv(u);
	float y = avgv(v);
	if (x < y)
		return 1;
	return -1;
}

void median_filter(Vector3 *src, Vector3 *dst,
	float *depth_buffer, float *roughness_buffer,
	Vector3 *normal_buffer, int frame_w, int frame_h,
	int region_x, int region_w, int region_h)
{
	for (int j = 0; j < region_h; j++)
		for (int i = region_x; i - region_x < region_w; i++) {

			int center_pixel_location = j * frame_w + i;

			Vector3 center_color = src[center_pixel_location];
			bool center_roughness = roughness_buffer[center_pixel_location];
/*
			if (center_roughness < 0.8) {
				filter_buffer[center_pixel_location] = color_buffer[center_pixel_location];
				continue;
			}
*/
			#define WINDOW_SIZE 5

			Vector3 samples[WINDOW_SIZE*WINDOW_SIZE];
			int num_samples = 0;

			for (int u = 0; u < WINDOW_SIZE; u++)
				for (int v = 0; v < WINDOW_SIZE; v++) {

					int g = j + u - WINDOW_SIZE / 2;
					int t = i + v - WINDOW_SIZE / 2;

					if (g < 0 || t < 0 || t >= frame_w || g >= frame_h)
						continue;

					int neighbor_pixel_location = g * frame_w + t;

					if (absf(depth_buffer[center_pixel_location] - depth_buffer[neighbor_pixel_location]) > 1)
						continue;

					bool neighbor_roughness = roughness_buffer[neighbor_pixel_location];
					if (absf(neighbor_roughness - center_roughness) > 0.2)
						continue;

					samples[num_samples++] = src[neighbor_pixel_location];
				}

			qsort(samples, num_samples, sizeof(Vector3), compare_vector_luminosity);

			dst[center_pixel_location] = samples[num_samples/2];
/*
			float error = (absf(center_color.x - result.x) + absf(center_color.y - result.y) + absf(center_color.z - result.z)) / 3;
			if (error > 0.05)
				filter_buffer[center_pixel_location] = result;
			else
				filter_buffer[center_pixel_location] = color_buffer[center_pixel_location];
*/
		}
	//memcpy(color_buffer, filter_buffer, sizeof(Vector3) * frame_w * frame_h);
}

void smooth_filter(Vector3 *src, Vector3 *dst,
	float *depth_buffer, float *roughness_buffer,
	Vector3 *normal_buffer, int frame_w, int frame_h,
	int region_x, int region_w, int region_h)
{
	for (int j = 0; j < region_h; j++)
		for (int i = region_x; i - region_x < region_w; i++) {

			int center_pixel_location = j * frame_w + i;

			Vector3 center_color = src[center_pixel_location];
			bool    center_roughness = roughness_buffer[center_pixel_location];

			float min_distance = 1;
			float max_distance = 8;
			float distance = clamp(depth_buffer[center_pixel_location], min_distance, max_distance) / max_distance;
			float sigma = 0.5 + (0.5 + 10 * center_roughness) / distance;

			// roughness=0, distance=0  -> sigma=0.5
			// roughness=1, distance=0  -> sigma=10
			// roughness=0, distance=15 -> sigma=0.5
			// roughness=1, distance=15 -> sigma=0.5

			//int min_winsize = 5;
			//int max_winsize = 15;
			//int winsize = (min_winsize + center_roughness * (max_winsize - min_winsize)) | 1;

			int winsize = 5;

			#define MAX_WINSIZE 32
			Vector3 samples[MAX_WINSIZE * MAX_WINSIZE];
			float   weights[MAX_WINSIZE * MAX_WINSIZE];
			int num_samples = 0;

			for (int u = 0; u < winsize; u++)
				for (int v = 0; v < winsize; v++) {

					int g = j + u - winsize / 2;
					int t = i + v - winsize / 2;

					if (g < 0 || t < 0 || t >= frame_w || g >= frame_h)
						continue;

					int neighbor_pixel_location = g * frame_w + t;

					if (absf(depth_buffer[center_pixel_location] - depth_buffer[neighbor_pixel_location]) > 1)
						continue;

					bool neighbor_roughness = roughness_buffer[neighbor_pixel_location];
					if (absf(neighbor_roughness - center_roughness) > 0.2)
						continue;

					Vector3 center_normal = normal_buffer[center_pixel_location];
					Vector3 neighbor_normal = normal_buffer[neighbor_pixel_location];
					float simil = dotv(center_normal, neighbor_normal);
					if (simil < 0.8)
						continue;

					float weight = exp(- (float) ((t - i)*(t - i) + (g - j)*(g - j)) / (2 * sigma * sigma));

					weights[num_samples] = weight;
					samples[num_samples] = src[neighbor_pixel_location];
					num_samples++;
				}

			int head = 0;
/*
			if (num_samples > 2) {
				qsort(samples, num_samples, sizeof(Vector3), compare_vector_luminosity);
				head++;
				num_samples -= 2;
			}
*/
			Vector3 result = {0, 0, 0};
			float weight_sum = 0;
			for (int i = head; i < num_samples; i++) {
				result.x += weights[i] * samples[i].x;
				result.y += weights[i] * samples[i].y;
				result.z += weights[i] * samples[i].z;
				weight_sum += weights[i];
			}

			result.x /= weight_sum;
			result.y /= weight_sum;
			result.z /= weight_sum;

			dst[center_pixel_location] = result;
		}
}

void clamp_filter(Vector3 *src, Vector3 *dst,
	float *depth_buffer, float *roughness_buffer,
	Vector3 *normal_buffer, int frame_w, int frame_h,
	int region_x, int region_w, int region_h)
{
	for (int j = 0; j < region_h; j++)
		for (int i = region_x; i - region_x < region_w; i++) {
			dst[j * frame_w + i].x = clamp(src[j * frame_w + i].x, 0, 1);
			dst[j * frame_w + i].y = clamp(src[j * frame_w + i].y, 0, 1);
			dst[j * frame_w + i].z = clamp(src[j * frame_w + i].z, 0, 1);

#if 1
			 // HDR tonemapping
			dst[j * frame_w + i].x = dst[j * frame_w + i].x / (dst[j * frame_w + i].x + 1);
			dst[j * frame_w + i].y = dst[j * frame_w + i].y / (dst[j * frame_w + i].y + 1);
			dst[j * frame_w + i].z = dst[j * frame_w + i].z / (dst[j * frame_w + i].z + 1);
#endif

#if 0
			// gamma correct
			dst[j * frame_w + i].x = pow(dst[j * frame_w + i].x, 1.0/2.2);
			dst[j * frame_w + i].y = pow(dst[j * frame_w + i].y, 1.0/2.2);
			dst[j * frame_w + i].z = pow(dst[j * frame_w + i].z, 1.0/2.2);
#endif
		}
}

uint64_t start_time_ns;
uint64_t start_time_cycles;

void update_frame(void)
{
	uint64_t frame_start = __rdtsc();

	if (frame_buffer_size_doesnt_match_window())
		realloc_frame_buffer();

	os_mutex_lock(&frame_mutex);

	completed = 0;
	global_frame_index++;
	for (int i = 0; i < num_columns; i++)
		os_condvar_signal(&start_work);

	while (completed < num_columns)
		os_condvar_wait(&completed_work, &frame_mutex, -1);
	
	os_mutex_unlock(&frame_mutex);

	color_buffer_frames++;

	bool smooth = false;
	if (smooth_index == smooth_limit) {
		smooth = true;
		smooth_index = 0;
		smooth_limit = (smooth_limit + 1) * SMOOTH_RATIO;
		printf("smoothing!\n");
	}
	smooth_index++;

	if (smooth) {
#if 1
		filter(color_buffer, filter_buffer_B, smooth_filter);
		memcpy(color_buffer, filter_buffer_B, frame_w * frame_h * sizeof(Vector3));
		filter(filter_buffer_B, filter_buffer_A, clamp_filter);
#else
		filter(color_buffer, filter_buffer_A, median_filter);
		filter(filter_buffer_A, filter_buffer_B, smooth_filter);
		memcpy(color_buffer, filter_buffer_B, frame_w * frame_h * sizeof(Vector3));
		filter(filter_buffer_B, filter_buffer_A, clamp_filter);
#endif
	} else {
		filter(color_buffer, filter_buffer_A, clamp_filter);
	}

	move_frame_to_the_gpu(frame_w, frame_h, filter_buffer_A);
	//move_frame_to_the_gpu(frame_w, frame_h, normal_buffer);

	uint64_t cycles_per_pixel = global_cycle_pixel_sum / global_pixel_count;

	uint64_t cycles_per_frame = __rdtsc() - frame_start;

	uint64_t current_time_cycles = __rdtsc();
	uint64_t current_time_ns = get_relative_time_ns();
	double cy2ns = (double) (current_time_ns - start_time_ns) / (current_time_cycles - start_time_cycles);

	printf("pixel -> %llu cycles (%f ns)\n", cycles_per_pixel, cycles_per_pixel * cy2ns);
	printf("frame -> %llu cycles (%f ns)\n", cycles_per_frame, cycles_per_frame * cy2ns);
}

void print_camera(void);

int main(int argc, char **argv)
{
	fprintf(stderr, "Started\n");

	char *scene_file;
	parse_arguments_or_exit(argc, argv, &num_columns, &scene_file);

	fprintf(stderr, "Parsed arguments\n");

	if (!parse_scene_file(scene_file, &scene)) {
		fprintf(stderr, "Couldn't parse scene\n");
		return -1;
	}

	fprintf(stderr, "Scene parsed\n");

	const char *faces[] = {
		[CF_RIGHT]  = "assets/skybox/right.jpg",
		[CF_LEFT]   = "assets/skybox/left.jpg",
		[CF_TOP]    = "assets/skybox/top.jpg",
		[CF_BOTTOM] = "assets/skybox/bottom.jpg",
		[CF_FRONT]  = "assets/skybox/front.jpg",
		[CF_BACK]   = "assets/skybox/back.jpg",
	};
	load_cubemap(&skybox, faces);

	fprintf(stderr, "Cubemap loaded\n");

#if 0
	{
		local_pixel_count = 0;
		local_pixel_cycles = 0;
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 1000; j++)
				for (int k = 0; k < 1000; k++) {
					float u = (float) j / 999;
					float v = (float) i / 999;
					Ray ray = ray_through_screen_at(u, v, 16.0f/9);
					//wyhash64_x = 0;
					pixel_inner(ray);
				}
		}
		printf("pixel -> %llu cycles\n", local_pixel_cycles / local_pixel_count);
		return 0;
	}
#endif

	startup_window_and_opengl_context_or_exit(256*4, 256*4, "Ray Tracing");

	fprintf(stderr, "Started windows and opengl context\n");

	start_ray_tracing_workers();
	start_filtering_workers();

	fprintf(stderr, "Workers started\n");

	start_time_ns = get_relative_time_ns();
	start_time_cycles = __rdtsc();

	for (bool exit = false; !exit; ) {

		//print_camera();

		for (;;) {

			double mouse_x;
			double mouse_y;
			int event = pop_event(&mouse_x, &mouse_y);
			if (event == EVENT_EMPTY) break;

			float speed = 0.5;
			switch (event) {
				case EVENT_CLOSE:
				case EVENT_PRESS_ESC:
				fprintf(stderr, "Exiting\n");
				exit = true;
				break;

				case EVENT_PRESS_W:
				case EVENT_AGAIN_W:
				move_camera(UP, speed);
				invalidate_accumulation();
				break;

				case EVENT_PRESS_A:
				case EVENT_AGAIN_A:
				move_camera(LEFT, speed);
				invalidate_accumulation();
				break;

				case EVENT_PRESS_S:
				case EVENT_AGAIN_S:
				move_camera(DOWN, speed);
				invalidate_accumulation();
				break;

				case EVENT_PRESS_D:
				case EVENT_AGAIN_D:
				move_camera(RIGHT, speed);
				invalidate_accumulation();
				break;

				case EVENT_MOVE_MOUSE:
				rotate_camera(mouse_x, mouse_y);
				invalidate_accumulation();
				break;

				case EVENT_PRESS_SPACE:
				screenshot();
				break;
			}
		}

		update_frame();
		draw_frame();
	}

	// Tell workers to stop evaluating frames
	invalidate_accumulation();

	stop_filtering_workers();
	stop_ray_tracing_workers();
	free_cubemap(&skybox);
	cleanup_window_and_opengl_context();
	return 0;
}

void parse_arguments_or_exit(int argc, char **argv, int *num_columns, char **scene_file)
{
	*scene_file = NULL;
	*num_columns = -1;
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--threads")) {
			i++;
			if (i == argc) {
				fprintf(stderr, "Error: --threads option is missing the count\n");
				exit(-1);
			}
			*num_columns = atoi(argv[i]);
			if (*num_columns == 0) {
				fprintf(stderr, "Error: Invalid count for --threads\n");
				exit(-1);
			}
		} else if (!strcmp(argv[i], "--scene")) {
			i++;
			if (i == argc) {
				fprintf(stderr, "Error: --scene option is missing the file path\n");
				exit(-1);
			}
			*scene_file = argv[i];
		} else {
			fprintf(stderr, "Warning: Ignoring option %s\n", argv[i]);
		}
	}
	if (*scene_file == NULL) {
		fprintf(stderr, "Error: No scene specified (you should use --scene <filename>)\n");
		exit(-1);
	}
	if (*num_columns < 0) {
		fprintf(stderr, "Error: Missing --threads <N> option\n");
		exit(-1);
	}
	if (*num_columns > MAX_COLUMNS)
		*num_columns = MAX_COLUMNS;
}

// Must be executed while holding the frame lock
void screenshot(void)
{
	// Choose a file name in the form "screenshot_X.png" where X
	// is an integer between 0 and 1000 that is not being used
	// already.
	char file[1<<12];
	int i = 0;
	while (i < 1000) {
		int k = snprintf(file, sizeof(file), "screenshot_%d.png", i);
		if (k < 0 || k >= (int) sizeof(file)) {
			fprintf(stderr, "Couldn't take screenshot (path buffer too small)\n");
			return;
		}
		FILE *stream = fopen(file, "rb");
		if (stream == NULL) {
			if (errno == ENOENT)
				break;
			fprintf(stderr, "Couldn't take screenshot (%s)\n", strerror(errno));
			return;
		}
		fclose(stream);
		i++;
	}

	// Convert the frame buffer from one float per pixel to one byte.
	uint8_t *converted = malloc(frame_w * frame_h * 3 * sizeof(uint8_t));
	if (converted == NULL) {
		fprintf(stderr, "Couldn't take screenshot (out of memory)\n");
	}
	for (int i = 0; i < frame_w * frame_h; i++) {
		converted[i * 3 + 0] = filter_buffer_A[i].x * 255;
		converted[i * 3 + 1] = filter_buffer_A[i].y * 255;
		converted[i * 3 + 2] = filter_buffer_A[i].z * 255;
	}

	stbi_flip_vertically_on_write(1);
	int ok = stbi_write_png(file, frame_w, frame_h, 3, converted, 0);

	free(converted);

	if (!ok)
		fprintf(stderr, "Could not take screenshot (write error)\n");
	else
		fprintf(stderr, "Took screenshot! (%s)\n", file);
}

/////////////////////////////////////////////////////////////////////////////
/// WORKER SYNCHRONIZATION                                                ///
/////////////////////////////////////////////////////////////////////////////

static bool workers_should_stop;
os_thread workers[MAX_COLUMNS];

bool quitting(void)
{
	return workers_should_stop;
}

void start_ray_tracing_workers(void)
{
	workers_should_stop = false;

	os_mutex_create(&frame_mutex);

	os_condvar_create(&start_work);
	os_condvar_create(&completed_work);

	for (int i = 0; i < num_columns; i++)
		os_thread_create(&workers[i], (void*) i, ray_tracing_worker);
}

void stop_ray_tracing_workers(void)
{
	os_mutex_lock(&frame_mutex);
	workers_should_stop = true;
	for (int i = 0; i < num_columns; i++)
		os_condvar_signal(&start_work);
	os_mutex_unlock(&frame_mutex);
	for (int i = 0; i < num_columns; i++)
		os_thread_join(workers[i]);

	os_condvar_delete(&start_work);
	os_condvar_delete(&completed_work);
}

///////////////////////////////////////////////////////////////////////////////////////

bool         filtering_workers_should_stop;
FilterFunc   filtering_func;
Vector3     *filtering_input;
Vector3     *filtering_output;
int          filtering_completed;
os_mutex_t   filtering_mutex;
os_condvar_t start_filtering;
os_condvar_t completed_filtering;
uint64_t     global_filter_index = 0;
os_thread    filtering_workers[MAX_COLUMNS];

os_threadreturn filtering_worker(void *arg);

void start_filtering_workers(void)
{
	filtering_workers_should_stop = false;

	os_mutex_create(&filtering_mutex);

	os_condvar_create(&start_filtering);
	os_condvar_create(&completed_filtering);

	for (int i = 0; i < num_columns; i++)
		os_thread_create(&filtering_workers[i], (void*) i, filtering_worker);
}

void stop_filtering_workers(void)
{
	os_mutex_lock(&filtering_mutex);
	filtering_workers_should_stop = true;
	for (int i = 0; i < num_columns; i++)
		os_condvar_signal(&start_filtering);
	os_mutex_unlock(&filtering_mutex);
	for (int i = 0; i < num_columns; i++)
		os_thread_join(workers[i]);

	os_condvar_delete(&start_work);
	os_condvar_delete(&completed_work);
}

void filter(Vector3 *input, Vector3 *output, FilterFunc func)
{
	filtering_func = func;
	filtering_input = input;
	filtering_output = output;

	os_mutex_lock(&filtering_mutex);

	filtering_completed = 0;
	global_filter_index++;
	for (int i = 0; i < num_columns; i++)
		os_condvar_signal(&start_filtering);

	while (filtering_completed < num_columns)
		os_condvar_wait(&completed_filtering, &filtering_mutex, -1);
	
	os_mutex_unlock(&filtering_mutex);
}

os_threadreturn filtering_worker(void *arg)
{
	int column_i = (int) arg;

	uint64_t local_filter_index = 0;

	os_mutex_lock(&filtering_mutex);
	for (;;) {

		while (local_filter_index == global_filter_index && !filtering_workers_should_stop)
			os_condvar_wait(&start_filtering, &filtering_mutex, -1);
		local_filter_index = global_filter_index;

		if (filtering_workers_should_stop) {
			filtering_completed++;
			os_condvar_signal(&completed_filtering);
			break;
		}

		os_mutex_unlock(&filtering_mutex);

		int column_w = frame_w / num_columns;
		int column_x = column_w * column_i;

		filtering_func(
			filtering_input,
			filtering_output,
			depth_buffer,
			roughness_buffer,
			normal_buffer, 
			frame_w,
			frame_h,
			column_x,
			column_w,
			frame_h);

		// Now we try publishing the changes
		os_mutex_lock(&filtering_mutex);
		filtering_completed++;
		os_condvar_signal(&completed_filtering);
	}
	os_mutex_unlock(&filtering_mutex);
}