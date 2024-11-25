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
#include <float.h>
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
	float f0;
	Vector3 emission;
} ExtraObjectData;

typedef void (*FilterFunc)(Vector3*, Vector3*, HitInfo*, int, int, int, int, int);

/////////////////////////////////////////////////////////////////////////////
/// GLOBAL VARIABLES                                                      ///
/////////////////////////////////////////////////////////////////////////////

#define MAX_COLUMNS 32
#define INIT_SMOOTH_LIMIT 1
#define SMOOTH_RATIO 1.01
#define FILTER 1
#define LIGHT_SAMPLING 1
#define IMPORTANCE_SAMPLING 1

// Parameters. These are set at startup and are
// considered constant after that.
int num_columns;

// The scene and background being rendered.
Scene   scene;
Cubemap skybox;

ExtraObjectData object_data[MAX_OBJECTS];

int color_buffer_frames = 0;

Vector3 *filter_buffer_A;
Vector3 *filter_buffer_B;

Vector3 *color_buffer;
HitInfo *hit_buffer;

// Size of the accumulation and frame buffers
int frame_w = 0;
int frame_h = 0;

// This guards the critical section around the accumulation buffer.
os_mutex_t frame_mutex;

int completed = 0;
os_condvar_t start_work;
os_condvar_t completed_work;

uint64_t global_frame_index = 0;

uint64_t global_cycle_pixel_sum = 0;
uint64_t global_pixel_count = 0;

uint64_t smooth_index = 0;
float    smooth_limit = INIT_SMOOTH_LIMIT;

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

Vector3 pixel(float x, float y, float aspect_ratio, HitInfo *first_hit);
void    update_frame(void);
void    invalidate_accumulation(void);

os_threadreturn worker(void *arg);

/////////////////////////////////////////////////////////////////////////////
/// IMPLEMENTATION                                                        ///
/////////////////////////////////////////////////////////////////////////////

void invalidate_accumulation(void)
{
	os_mutex_lock(&frame_mutex);
	memset(color_buffer, 0, sizeof(Vector3) * frame_w * frame_h);
	memset(hit_buffer, 0, sizeof(HitInfo) * frame_w * frame_h);
	color_buffer_frames = 0;
	smooth_index = 0;
	smooth_limit = INIT_SMOOTH_LIMIT;
	os_mutex_unlock(&frame_mutex);
}

static Vector3 uniform_sample(Vector3 n)
{
	float u1 = random_float();
	float u2 = random_float();

	float r = sqrt(1.0 - u1 * u1);
	float phi = 2.0 * M_PI * u2;

	Vector3 B = normalize(cross(n, (Vector3) {0, 1, 1}));
	Vector3 T = cross(B, n);

	Vector3 R;
	R.x = r * sin(phi) * B.x + u1 * n.x + r * cos(phi) * T.x;
	R.y = r * sin(phi) * B.y + u1 * n.y + r * cos(phi) * T.y;
	R.z = r * sin(phi) * B.z + u1 * n.z + r * cos(phi) * T.z;
	return normalize(R);
}

static Vector3 cosine_sample(Vector3 n, float max_radius)
{
	float u1 = random_float();
	float u2 = random_float();

	float r = sqrt(u1) * max_radius;
	float theta = 2.0 * M_PI * u2;

	Vector3 B = normalize(cross(n, (Vector3) {0, 1, 1}));
	Vector3 T = cross(B, n);

	Vector3 R;
	R.x = r * sin(theta) * B.x + sqrt(1.0 - u1) * n.x + r * cos(theta) * T.x;
	R.y = r * sin(theta) * B.y + sqrt(1.0 - u1) * n.y + r * cos(theta) * T.y;
	R.z = r * sin(theta) * B.z + sqrt(1.0 - u1) * n.z + r * cos(theta) * T.z;
	return normalize(R);
}

// Generate a random number in [0, 1] and return its probability.
//
// This is the pdf:
//
//   |        /\ 
//   |       /: \ 
//   |      / :  \ 
//   |_____/  :   \_____ _ _ _
//         0  :  1
//            :
//           peak
//
// It's only symmetric when peak is 0.5
float random_with_peak(float avg, float *p)
{
	avg = clamp(avg, 0.0001, 1 - 0.0001);
	float u = random_float();
	if (u < avg)
		u = sqrt(u) * avg / sqrt(avg);
	else
		u = 1 - sqrt(u - avg) * (1 - avg) / sqrt(1 - avg);

	if (u < avg)
		*p = 2 * u / avg;
	else
		*p = 2 * (1 - u) / (1 - avg);

	return u;
}

// Generate a random number in [0, 1] and return its probability.
//
// This is the pdf:
//
//   |        /\ 
//   |       /: \ 
//   |      / :  |
//   |_____/  :  |_____ _ _ _
//         0  :  1
//            :
//           peak
//
// Symmetric around "peak" but truncated at 1
float random_with_peak_trunc(float peak, float *p)
{
	peak = 1 - peak;
	if (peak < 0.5 || peak > 1) {
		printf("invalid peak %f\n", peak);
		abort();
	}
	peak = clamp(peak, 0.5 + 0.0001, 1 - 0.0001);
	float u;
	do {
		u = random_float();
		if (u < 0.5)
			u = sqrt(u) * 0.5 / sqrt(0.5);
		else
			u = 1 - sqrt(u - 0.5) * 0.5 / sqrt(0.5);
		u *= 2 * peak;
	} while (u > 1);

	if (u < peak)
		*p = 0.5 * u / ((peak - 0.5) * (peak - 0.5));
	else
		*p = (peak - 0.5 * u) / ((peak - 0.5) * (peak - 0.5));

	u = 1 - u;
	return u;
}

// Generate a random direction d such that dot(d,n)>0. Results
// are biased towards the "b" direction. The probability of the
// generated vector is returned through "p".
static Vector3 sample_normal_hemisphere_with_bias(Vector3 n, Vector3 b, float *p)
{
	n = normalize(n);
	b = normalize(b);

	float phi;
	float theta;
	float p1;
	float p2;

	Vector3 s;

	float NoB = dotv(n, b);
	if (NoB < 0) {
		printf("We are assuming n-b point in the same general direction\n");
		abort();
	}

	// Vector orthogonal to the normal
	s.x = b.x - n.x * NoB;
	s.y = b.y - n.y * NoB;
	s.z = b.z - n.z * NoB;
	s = normalize(s);

	// Calculate the angle on the n-s plane
	float phi_bias = M_PI/2 - acos(NoB);
	phi = random_with_peak_trunc(phi_bias / M_PI, &p1) * M_PI;

	// Angle on the n-t plane
	theta = random_with_peak(0.5, &p2) * M_PI;

	// Now calculate the third vector for an R3 basis
	Vector3 t = cross(n, s);
	t = normalize(t);

	// Coordinates in the (n,s,t) base
	float w0 = sin(phi);
	float w2 = cos(phi) * cos(theta);
	float w1 = cos(phi) * sin(theta);

	Vector3 r;
	r.x = n.x * w0 + s.x * w1 + t.x * w2;
	r.y = n.y * w0 + s.y * w1 + t.y * w2;
	r.z = n.z * w0 + s.z * w1 + t.z * w2;
	r = normalize(r);

	*p = p1 * p2;
	return r;
}

Ray move_ray_origin_forward(Ray ray, float delta_t)
{
	ray.origin.x += ray.direction.x * delta_t;
	ray.origin.y += ray.direction.y * delta_t;
	ray.origin.z += ray.direction.z * delta_t;
	return ray;
}

float update_average(float avg, float sample, int num)
{
	return (avg * num + sample) / (num + 1);
}

Vector3 update_average_vec(Vector3 avg, Vector3 sample, int num)
{
	Vector3 R;
	R.x = update_average(avg.x, sample.x, num);
	R.y = update_average(avg.y, sample.y, num);
	R.z = update_average(avg.z, sample.z, num);
	return R;
}

Vector3 pixel(float x, float y, float aspect_ratio, HitInfo *first_hit)
{
	Ray ray = ray_through_screen_at(x, y, aspect_ratio);

	int light_index = -1;
	Vector3 light_origin;
	for (int i = 0; i < scene.num_objects; i++) {
		if (scene.objects[i].material.emission_power > 0) {
			light_index = i;
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
	int bounces = 5;

	HitInfo default_hit;
	default_hit.normal = scalev(ray.direction, -1);
	default_hit.distance = FLT_MAX;
	default_hit.object = -1;

	*first_hit = default_hit;

	for (int i = 0; i < bounces; i++) {

		// Find the next collision
		HitInfo hit = trace_ray(ray, &scene);
		if (hit.object == -1) {
			// The ray flew straight out of the scene!
			//
			// The sky is sampled here. You can change the sky color here
			// if you want:
			//     Vector3 sky_color = {0.6, 0.7, 0.9};
			//     Vector3 sky_color = {0, 0, 0};
			//     Vector3 sky_color = {1, 1, 1};
			Vector3 sky_color = sample_cubemap(&skybox, normalize(ray.direction));
			result = combine(result, mulv(sky_color, contrib), 1, 1);
			break;
		}
		if (i == 0) *first_hit = hit;

		Material material = scene.objects[hit.object].material;

		Vector3 v = normalize(ray.direction);
		Vector3 n = normalize(hit.normal);
		Vector3 o = hit.point;
		float   f0 = object_data[hit.object].f0;
		Vector3 em = object_data[hit.object].emission;

		result.x += contrib.x * em.x;
		result.y += contrib.y * em.y;
		result.z += contrib.z * em.z;

		float NoV = clamp(-dotv(n, v), 0, 1);

		float F = f0 + (1 - f0) * pow(1 - NoV, 5);

		Vector3 shade;
		Vector3 out_dir;
		if (random_float() < F) {

			// Specular ray

			Vector3 r = reflect(v, n);
			out_dir = cosine_sample(r, material.roughness);
			if (dotv(out_dir, n) < 0) out_dir = scalev(out_dir, -1);

			shade = material.albedo;

		} else {

			// Diffuse ray

			float factor;

#if LIGHT_SAMPLING
			{
				if (light_index < 0) {

					// No lights in the scene. Fallback on cosine sampling
					out_dir = cosine_sample(n, 1);
					factor = 1 / clamp(dotv(out_dir, n), FLT_MIN, 1);

				} else {

					Vector3 l;
					l.x = light_origin.x - o.x;
					l.y = light_origin.y - o.y;
					l.z = light_origin.z - o.z;
					l = normalize(l);

					// When the light is not in the general direction of the normal,
					// we can't do light sampling. To avoid abrupt transitions we mix
					// cosine and light sampling.
					if (random_float() < dotv(l, n)) {
						float p;
						out_dir = sample_normal_hemisphere_with_bias(n, l, &p);
						factor =  1 / clamp(p, 0.0001, 1);
					} else {
						out_dir = cosine_sample(n, 1);
						factor = 1 / clamp(dotv(out_dir, n), 0.0001, 1);
					}
				}
			}
#elif IMPORTANCE_SAMPLING
			{
				out_dir = cosine_sample(n, 1);
				factor = 1 / clamp(dotv(out_dir, n), FLT_MIN, 1);
			}
#else
			{
				out_dir = uniform_sample(n);
				factor = 2;
			}
#endif

			factor *= clamp(dotv(out_dir, n), 0, 1);
			shade = scalev(material.albedo, factor);
		}

		contrib.x *= shade.x;
		contrib.y *= shade.y;
		contrib.z *= shade.z;

		ray.origin = o;
		ray.direction = out_dir;
		ray = move_ray_origin_forward(ray, 0.001);
	}

	// Avoid fireflies
	result.x = clamp(result.x, 0, 10);
	result.y = clamp(result.y, 0, 10);
	result.z = clamp(result.z, 0, 10);

	return result;
}

os_threadreturn ray_tracing_worker(void *arg)
{
	Vector3 *local_color_buffer = NULL;
	HitInfo *local_hit_buffer = NULL;

	int column_i = (int) arg;

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
		int column_w = frame_w / num_columns;
		cached_frame_w = frame_w;
		cached_frame_h = frame_h;

		os_mutex_unlock(&frame_mutex);

		// We need to resize
		if (resize) {

			free(local_color_buffer);
			local_color_buffer = malloc(sizeof(Vector3) * column_w * cached_frame_h);
			if (!local_color_buffer) abort();

			free(local_hit_buffer);
			local_hit_buffer = malloc(sizeof(HitInfo) * column_w * cached_frame_h);
			if (!local_hit_buffer) abort();
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

				// Skip background pixels that were already calculated
				if (color_buffer_frames > 0 && local_hit_buffer[j * column_w + i].object == -1)
					continue;

				uint64_t start_time = __rdtsc();
				local_color_buffer[j * column_w + i] = pixel(u, v, aspect_ratio, &local_hit_buffer[j * column_w + i]);
				local_cycle_pixel_sum += __rdtsc() - start_time;
				local_pixel_count++;

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

				color_buffer[dst_index] = update_average_vec(color_buffer[dst_index], local_color_buffer[src_index], color_buffer_frames);

				hit_buffer[dst_index] = local_hit_buffer[src_index];
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

	free(hit_buffer);
	hit_buffer = malloc(sizeof(HitInfo) * frame_w * frame_h);
	if (!hit_buffer) {
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
	memset(hit_buffer,  0, sizeof(HitInfo) * frame_w * frame_h);
}

bool frame_buffer_size_doesnt_match_window(void)
{
	return frame_w != get_screen_w() || frame_h != get_screen_h();
}

int compare_colors(const void *a, const void *b)
{
	Vector3 *a2 = (Vector3*) a;
	Vector3 *b2 = (Vector3*) b;
	float asum = a2->x + a2->y + a2->z;
	float bsum = b2->x + b2->y + b2->z;
	if (asum > bsum)
		return 1;
	return -1;
}

void smooth_filter(Vector3 *src, Vector3 *dst, HitInfo *hit_buffer,
	int frame_w, int frame_h, int region_x, int region_w, int region_h)
{
	for (int j = 0; j < region_h; j++)
		for (int i = region_x; i - region_x < region_w; i++) {

			int center_pixel_location = j * frame_w + i;

			Vector3 center_color = src[center_pixel_location];
			HitInfo center_hit = hit_buffer[center_pixel_location];
			if (center_hit.object < 0) {
				dst[center_pixel_location] = src[center_pixel_location];
				continue; // Skip background pixels
			}
			Material material = scene.objects[center_hit.object].material;

			bool uniform_weights;
			float sigma;
			int winsize;

			winsize = 3;
			sigma = 1;

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

					HitInfo neighbor_hit = hit_buffer[neighbor_pixel_location];

					if (neighbor_hit.object != center_hit.object || dotv(center_hit.normal, neighbor_hit.normal) < 0.8)
						continue;

					float weight = exp(- (float) ((t - i)*(t - i) + (g - j)*(g - j)) / (2 * sigma * sigma));

					Vector3 sample = src[neighbor_pixel_location];

					weights[num_samples] = weight;
					samples[num_samples] = sample;
					num_samples++;
				}

			if (num_samples == 0)
				continue;

			int head = 0;
			int skip = 0;
			if (skip > 0) {
				qsort(samples, num_samples, sizeof(Vector3), compare_colors);
				if (num_samples > 2 * skip) {
					head = skip;
					num_samples -= 2 * skip;
				} else {
					head = num_samples/2;
					num_samples = 1;
				}
			}

			Vector3 result = {0, 0, 0};
			float weight_sum = 0;
			for (int i = head; i - head < num_samples; i++) {
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

Vector3 white;

float map_tone_channel(float C, float W)
{
	return C * (1 + C / (W * W)) / (1 + C);
}

Vector3 map_tone(Vector3 C, Vector3 W)
{
	Vector3 R;
	R.x = map_tone_channel(C.x, W.x);
	R.y = map_tone_channel(C.y, W.y);
	R.z = map_tone_channel(C.z, W.z);
	return R;
}

void clamp_filter(Vector3 *src, Vector3 *dst, HitInfo *hit_buffer,
	int frame_w, int frame_h, int region_x, int region_w, int region_h)
{
	for (int j = 0; j < region_h; j++)
		for (int i = region_x; i - region_x < region_w; i++) {
			dst[j * frame_w + i] = map_tone(src[j * frame_w + i], white);
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

#if FILTER
	bool smooth = false;
	if (smooth_index >= smooth_limit) {
		smooth = true;
		smooth_index = 0;
		smooth_limit *= SMOOTH_RATIO;
		printf("smoothing! (limit %f)\n", smooth_limit);
	} else {
		smooth_index++;
	}
#else
	bool smooth = false;
#endif

	white.x = 0;
	white.y = 0;
	white.z = 0;
	for (int j = 0; j < frame_h; j++)
		for (int i = 0; i < frame_w; i++) {
			Vector3 color = color_buffer[j * frame_w + i];
			white.x = maxf(white.x, color.x);
			white.y = maxf(white.y, color.y);
			white.z = maxf(white.z, color.z);
		}

	if (smooth) {
		filter(color_buffer, filter_buffer_B, smooth_filter);
		memcpy(color_buffer, filter_buffer_B, frame_w * frame_h * sizeof(Vector3));
		filter(filter_buffer_B, filter_buffer_A, clamp_filter);
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

	for (int i = 0; i < scene.num_objects; i++) {
		Material m = scene.objects[i].material;
		object_data[i].emission = m.emission_color;
		object_data[i].emission.x *= m.emission_power;
		object_data[i].emission.y *= m.emission_power;
		object_data[i].emission.z *= m.emission_power;
		object_data[i].f0 = 0.16 * m.reflectance * m.reflectance * (1 - m.metallic) + avgv(m.albedo) * m.metallic;
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
			hit_buffer,
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