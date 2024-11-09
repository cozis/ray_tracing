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

/////////////////////////////////////////////////////////////////////////////
/// GLOBAL VARIABLES                                                      ///
/////////////////////////////////////////////////////////////////////////////

#define MAX_COLUMNS 32

// Parameters. These are set at startup and are
// considered constant after that.
int num_columns;
int init_scale;

// The scene and background being rendered.
Scene   scene;
Cubemap skybox;

// Any time the accumulation buffer is reset or
// resized, this is incremented.
_Atomic uint32_t accum_generation = 0;

// This is the "accumulation buffer". Workers evaluate
// pixel colors in parallel and sum their results in here.
// When the main thread needs to draw a new frame it takes
// these values and divides them by the frame count, averaging
// the results of multiple frames.
Vector3 *accum = NULL;

// This is the "frame buffer". It's only accessed by the
// main buffer to store the averaged values of the accumulation
// buffer before sending them to the GPU.
Vector3 *frame = NULL;

// Size of the accumulation and frame buffers
int frame_w = 0;
int frame_h = 0;

// This guards the critical section around the accumulation buffer.
os_mutex_t frame_mutex;

// One condition variable per column. Any time new information
// is added to the accumulation buffer the condition of the
// associated column is signaled
os_condvar_t accum_conds[MAX_COLUMNS];

// Counters that indicate how much information each column
// is storing. An integer value of N means N full frames have 
// been accumulated. Lower resolution frames contribute lower
// values (half resolution weighs 0.25).
float accum_counts[MAX_COLUMNS];

/////////////////////////////////////////////////////////////////////////////
/// FUNCTION PROTOTYPES                                                   ///
/////////////////////////////////////////////////////////////////////////////

void    start_workers(void);
void    stop_workers(void);

bool    quitting(void);
void    screenshot(void);
void    parse_arguments_or_exit(int argc, char **argv, int *num_columns, int *init_scale, char **scene_file);

Vector3 pixel(float x, float y, float aspect_ratio);
void    update_frame(void);
float   render_column(Vector3 *data, int scale, int column_w, int column_i, int frame_w, int frame_h, uint64_t cached_generation);
void    invalidate_accumulation(void);

os_threadreturn worker(void *arg);

/////////////////////////////////////////////////////////////////////////////
/// IMPLEMENTATION                                                        ///
/////////////////////////////////////////////////////////////////////////////

// Resets the current frame and accumulation buffers and tells
// every worker to drop what they are doing and start again.
void invalidate_accumulation(void)
{
	os_mutex_lock(&frame_mutex);
	for (int i = 0; i < num_columns; i++)
		accum_counts[i] = 0;
	atomic_fetch_add(&accum_generation, 1);
	memset(accum, 0, sizeof(Vector3) * frame_w * frame_h);
	memset(frame, 0, sizeof(Vector3) * frame_w * frame_h);
	os_mutex_unlock(&frame_mutex);
}

Vector3 fresnel_schlick(float u, Vector3 f0)
{
	return combine(f0, combine(vec_from_scalar(1.0), f0, 1, -1), 1, pow(1.0 - u, 5.0));
}

static _Atomic uint64_t pixel_cycles = 0;
static _Atomic uint64_t pixel_count  = 0;

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

Vector3 pixel_inner(Ray in_ray)
{
	uint64_t start_time = __rdtsc();

	// Find a light source. This is kind of lazy as we should
	// sample every light source in the scene.
	int light_index = -1;
	float light_sample_weight = 0.05;
	Vector3 weighted_light_emission = {0, 0, 0};
	for (int i = 0; i < scene.num_objects; i++) {
		Material material = scene.objects[i].material;
		if (material.emission_power > 0) {
			light_index = i;
			weighted_light_emission.x += material.emission_color.x * material.emission_power * light_sample_weight;
			weighted_light_emission.y += material.emission_color.y * material.emission_power * light_sample_weight;
			weighted_light_emission.z += material.emission_color.z * material.emission_power * light_sample_weight;
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

		// Sample the light source
		//
		// Because we are only calculating on ray per pixel each frame, the impact
		// if light sources is greatly underestimated. In this loop we try hitting
		// light explicitly.
		bool light_sampled = false;
		if (light_index != -1) {

			Vector3 hp = hit.point;
			Vector3 hn = hit.normal;

			Object *object = &scene.objects[light_index];
			Vector3 origin = origin_of(*object);

			// Direction from the current collusion point to the light source
			Vector3 dir_to_light;
			dir_to_light.x = origin.x - hp.x;
			dir_to_light.y = origin.y - hp.y;
			dir_to_light.z = origin.z - hp.z;

			float spread = 0.5;

			Vector3 rand_dir;
			rand_dir.x = 2 * (float) wyhash64() / UINT64_MAX - 1;
			rand_dir.y = 2 * (float) wyhash64() / UINT64_MAX - 1;
			rand_dir.z = 2 * (float) wyhash64() / UINT64_MAX - 1;

			float dot = rand_dir.x * hn.x + rand_dir.y * hn.y + rand_dir.z * hn.z;
			if (dot < 0) {
				rand_dir.x = -rand_dir.x;
				rand_dir.y = -rand_dir.y;
				rand_dir.z = -rand_dir.z;
			}

			Vector3 sample_dir;
			sample_dir.x = spread * rand_dir.x + dir_to_light.x;
			sample_dir.y = spread * rand_dir.y + dir_to_light.y;
			sample_dir.z = spread * rand_dir.z + dir_to_light.z;

			Ray sample_ray;
			float eps = 0.001;
			sample_ray.direction.x = sample_dir.x;
			sample_ray.direction.y = sample_dir.y;
			sample_ray.direction.z = sample_dir.z;
			sample_ray.origin.x = hp.x + eps * sample_dir.x;
			sample_ray.origin.y = hp.y + eps * sample_dir.y;
			sample_ray.origin.z = hp.z + eps * sample_dir.z;

			HitInfo hit2 = trace_ray(sample_ray, &scene);
			if (hit2.object == light_index)
				light_sampled = true;
		}

		Material material = scene.objects[hit.object].material;

		Vector3 v;
		v.x = -in_ray.direction.x;
		v.y = -in_ray.direction.y;
		v.z = -in_ray.direction.z;

		Vector3 n = hit.normal;

		float NoV = n.x * v.x + n.y * v.y + n.z * v.z;
		if (NoV < 0)
			NoV = 0;
		else {
			if (NoV > 1)
				NoV = 1;
		}

		float f0 = 0.16 * material.reflectance * material.reflectance * (1 - material.metallic) + avgv(material.albedo) * material.metallic;
		float F = f0 + (1 - f0) * pow(1 - NoV, 5);

		Vector3 rand_dir;
		{
			rand_dir.x = 2 * (float) wyhash64() / UINT64_MAX - 1;
			rand_dir.y = 2 * (float) wyhash64() / UINT64_MAX - 1;
			rand_dir.z = 2 * (float) wyhash64() / UINT64_MAX - 1;
			if (rand_dir.x * n.x + rand_dir.y * n.y + rand_dir.z * n.z < 0) {
				rand_dir.x = -rand_dir.x;
				rand_dir.y = -rand_dir.y;
				rand_dir.z = -rand_dir.z;
			}
		}

		result.x += contrib.x * material.emission_color.x * material.emission_power;
		result.y += contrib.y * material.emission_color.y * material.emission_power;
		result.z += contrib.z * material.emission_color.z * material.emission_power;

		Vector3 out_dir;
		if (material.metallic > 0.001 || random_float() <= F) {

			// Specular ray

			Vector3 reflect_dir;
			float tmp = dotv(n, v);
			reflect_dir.x = n.x * 2 * tmp - v.x;
			reflect_dir.y = n.y * 2 * tmp - v.y;
			reflect_dir.z = n.z * 2 * tmp - v.z;

			out_dir.x = rand_dir.x * material.roughness + reflect_dir.x;
			out_dir.y = rand_dir.y * material.roughness + reflect_dir.y;
			out_dir.z = rand_dir.z * material.roughness + reflect_dir.z;

		} else {

			// Diffuse ray
			out_dir = rand_dir;

			contrib.x *= material.albedo.x * (1 - material.metallic);
			contrib.y *= material.albedo.y * (1 - material.metallic);
			contrib.z *= material.albedo.z * (1 - material.metallic);
		}

		Ray out_ray;
		out_ray.direction.x = out_dir.x;
		out_ray.direction.y = out_dir.y;
		out_ray.direction.z = out_dir.z;
		out_ray.origin.x = hit.point.x + out_dir.x * 0.001;
		out_ray.origin.y = hit.point.y + out_dir.y * 0.001;
		out_ray.origin.z = hit.point.z + out_dir.z * 0.001;

		if (light_sampled) {

			result.x += contrib.x * weighted_light_emission.x;
			result.y += contrib.y * weighted_light_emission.y;
			result.z += contrib.z * weighted_light_emission.z;

			contrib.x *= 1 - light_sample_weight;
			contrib.y *= 1 - light_sample_weight;
			contrib.z *= 1 - light_sample_weight;
		}

		in_ray = out_ray;
	}

	// Saturate the result so it's a valid color
	result.x = clamp(result.x, 0, 1);
	result.y = clamp(result.y, 0, 1);
	result.z = clamp(result.z, 0, 1);

	uint64_t end_time = __rdtsc() - start_time;
	atomic_fetch_add(&pixel_cycles, end_time);
	atomic_fetch_add(&pixel_count, 1);

	return result;
}

Vector3 pixel(float x, float y, float aspect_ratio)
{
	assert(!isnan(aspect_ratio));

	Ray in_ray = ray_through_screen_at(x, y, aspect_ratio);
	assert(!isnanv(in_ray.direction));

	return pixel_inner(in_ray);
}

float render_column(Vector3 *data, int scale, int column_w, int column_i, int frame_w, int frame_h, uint64_t cached_generation)
{
	// Since we're rendering at lower resolution, the weight of the
	// pixels we produce is also reduced.
	float scale2inv = 1.0f / (scale * scale);

	int column_x = column_w * column_i;
	float aspect_ratio = (float) frame_w / frame_h;

	// Just lower resolution version of each variable 
	int lowres_frame_w = frame_w / scale;
	int lowres_frame_h = frame_h / scale;
	int lowres_column_w = column_w / scale + 1;
	int lowres_column_x = column_x / scale;

	// Iterate over each low resolution pixel
	for (int j = 0; j < lowres_frame_h; j++) {
		for (int i = 0; i < lowres_column_w; i++) {

			float u = (float) (lowres_column_x + i) / (lowres_frame_w - 1);
			float v = (float) j  / (lowres_frame_h - 1);
			u = 1 - u;
			v = 1 - v;

			// Now copy the value of the single low resolution
			// pixel into a square of high resolution pixels
			int tile_w = scale;
			int tile_h = scale;
			if (tile_w > column_w - i * scale)
				tile_w = column_w - i * scale;

			Vector3 color = pixel(u, v, aspect_ratio);

			for (int g = 0; g < tile_h; g++)
				for (int t = 0; t < tile_w; t++) {
					int pixel_index = (j * scale + g) * column_w + (i * scale + t);
					assert(pixel_index >= 0 && pixel_index < column_w * frame_h);
					data[pixel_index] = scalev(color, 1);
				}
		}
		// We are done calculating a row of pixels!
		
		// If the frame has been invalidated we need to
		// exit and try again as soon as possible
		if (cached_generation != atomic_load(&accum_generation))
			break;
	}

	// Return the weight of the current column
	return scale2inv;
}

os_threadreturn worker(void *arg)
{
	// How many information is contained in the column buffer
	float column_data_weight = 0;

	// The actual pixels
	Vector3 *column_data = NULL;

	// The screen is divided in "num_columns" columns
	int column_i = (int) arg;
	int column_w;

	// Workers need to know the frame size while evaluating pixel
	// values. Since the frame size may change at any time, threads
	// cache their value.
	int cached_frame_w;
	int cached_frame_h;

	// Generation counter of the frame buffer when the worker
	// started producing a new frame. If the camera moves in the
	// or something else causing the frame buffer to be reset, this
	// will let the worker know the information needs to be thrown
	// away. 
	uint64_t cached_generation;

	// This value determines the resolution at which pixels are
	// evaluated. For scale=1 the image is full size. For scale=2
	// the image size is halved (along both axis). When a worker
	// evaluates a frame it starts at the lowest resolution "init_scale"
	// and after each succesfull paint it doubles the resolution
	int scale = init_scale;

	os_mutex_lock(&frame_mutex);
	while (!quitting()) {

		// Cache data and check if we need to resize the column buffer
		bool resize = false;
		if (column_data == NULL || cached_generation != atomic_load(&accum_generation))
			resize = true;
		column_w = frame_w / num_columns;
		cached_frame_w = frame_w;
		cached_frame_h = frame_h;
		cached_generation = atomic_load(&accum_generation);
		os_mutex_unlock(&frame_mutex);

		// We need to resize
		if (resize) {
			free(column_data);
			column_data = malloc(sizeof(Vector3) * column_w * cached_frame_h);
			if (!column_data) abort();
		}

		// Trace rays for each pixel in the column
		column_data_weight += render_column(column_data, scale, column_w, column_i, cached_frame_w, cached_frame_h, cached_generation);

		// Now we try publishing the changes
		os_mutex_lock(&frame_mutex);

		if (cached_generation == atomic_load(&accum_generation)) {
			// Frame didn't change its size while we were evaluating the column

			// This loop basically copies the pixel colors from the column buffer to
			// the frame buffer.
			for (int j = 0; j < frame_h; j++)
				for (int i = 0; i < column_w; i++) {	
					int column_x = column_w * column_i;
					int src_index = j * column_w + i;
					int dst_index = j * frame_w + (i + column_x);
					assert(src_index >= 0 && src_index < column_w * cached_frame_h);
					assert(dst_index >= 0 && dst_index < cached_frame_w * cached_frame_h);
					accum[dst_index] = combine(accum[dst_index], column_data[src_index], 1, 1.0f / (scale * scale));
				}
			accum_counts[column_i] += column_data_weight;

			// Let the main thread know there are new pixels
			os_condvar_signal(&accum_conds[column_i]);

			// We painted succesfully so we can render at double the resolution next time
			if (scale > 1)
				scale >>= 1;

		} else {
			// Data was invalidated. We need to go back and render at low res
			scale = init_scale;
		}

		// Either way we need to reset the column data now
		column_data_weight = 0;
	}
	os_mutex_unlock(&frame_mutex);
}

void realloc_frame_buffer(void)
{
	frame_w = get_screen_w();
	frame_h = get_screen_h();

	if (frame) free(frame);
	if (accum) free(accum);

	frame = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!frame) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	accum = malloc(sizeof(Vector3) * frame_w * frame_h);
	if (!accum) {
		printf("OUT OF MEMORY\n");
		abort();
	}

	for (int i = 0; i < num_columns; i++)
		accum_counts[i] = 0;
		
	memset(accum, 0, sizeof(Vector3) * frame_w * frame_h);
	memset(frame, 0, sizeof(Vector3) * frame_w * frame_h);

	atomic_fetch_add(&accum_generation, 1);
}

bool frame_buffer_size_doesnt_match_window(void)
{
	return frame_w != get_screen_w() || frame_h != get_screen_h();
}

void update_frame(void)
{
	os_mutex_lock(&frame_mutex);

	if (frame_buffer_size_doesnt_match_window())
		realloc_frame_buffer();

	int column_w = frame_w / num_columns;

	// Wait for the workers to produce a frame
	// (each worker produces a column)
	for (int i = 0; i < num_columns; i++) {
		while (accum_counts[i] < 0.0001)
			os_condvar_wait(&accum_conds[i], &frame_mutex, -1);
	}

	// Copy pixels from the accumulation buffer to the frame buffer
	for (int j = 0; j < frame_h; j++)
		for (int i = 0; i < frame_w; i++) {

			float u = (float) i / (frame_w - 1);
			float v = (float) j / (frame_h - 1);
			u = 1 - u;
			v = 1 - v;

			int pixel_index = j * frame_w + i;
			frame[pixel_index] = scalev(accum[pixel_index], 1.0f / accum_counts[i / column_w]);
		}

	move_frame_to_the_gpu(frame_w, frame_h, frame);

	os_mutex_unlock(&frame_mutex);
}

int main(int argc, char **argv)
{
/*
	{
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
		uint64_t pixel_count_2  = atomic_load(&pixel_count);
		uint64_t pixel_cycles_2 = atomic_load(&pixel_cycles);
		printf("pixel -> %llu cycles\n", pixel_cycles_2 / pixel_count_2);
		return 0;
	}
*/

	fprintf(stderr, "Started\n");

	char *scene_file;
	parse_arguments_or_exit(argc, argv, &num_columns, &init_scale, &scene_file);

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

	startup_window_and_opengl_context_or_exit(2 * 640, 2 * 480, "Ray Tracing");

	fprintf(stderr, "Started windows and opengl context\n");

	start_workers();

	fprintf(stderr, "Workers started\n");

	for (bool exit = false; !exit; ) {

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

		uint64_t pixel_count_2  = atomic_load(&pixel_count);
		uint64_t pixel_cycles_2 = atomic_load(&pixel_cycles);
		printf("pixel -> %llu cycles\n", pixel_cycles_2 / pixel_count_2);
	}

	// Tell workers to stop evaluating frames
	invalidate_accumulation();

	stop_workers();
	free_cubemap(&skybox);
	cleanup_window_and_opengl_context();
	return 0;
}

void parse_arguments_or_exit(int argc, char **argv, int *num_columns, int *init_scale, char **scene_file)
{
	*scene_file = NULL;
	*num_columns = -1;
	*init_scale = 8;
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--init-scale")) {
			i++;
			if (i == argc) {
				fprintf(stderr, "Error: --threads option is missing the count\n");
				exit(-1);
			}
			*init_scale = atoi(argv[i]);
			if (*init_scale != 1 && *init_scale != 2 && *init_scale != 4 && *init_scale != 8 && *init_scale != 16) {
				fprintf(stderr, "Error: Invalid value for --init-scale. It must be a power of 2 between 1 and 16 (included)\n");
				exit(-1);
			}
		} else if (!strcmp(argv[i], "--threads")) {
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
		converted[i * 3 + 0] = frame[i].x * 255;
		converted[i * 3 + 1] = frame[i].y * 255;
		converted[i * 3 + 2] = frame[i].z * 255;
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

void start_workers(void)
{
	workers_should_stop = false;

	os_mutex_create(&frame_mutex);

	for (int i = 0; i < num_columns; i++)
		os_condvar_create(&accum_conds[i]);

	for (int i = 0; i < num_columns; i++)
		os_thread_create(&workers[i], (void*) i, worker);
}

void stop_workers(void)
{
	os_mutex_lock(&frame_mutex);
	workers_should_stop = true;
	os_mutex_unlock(&frame_mutex);
	for (int i = 0; i < num_columns; i++)
		os_thread_join(workers[i]);

	for (int i = 0; i < num_columns; i++)
		os_condvar_delete(&accum_conds[i]);
}
