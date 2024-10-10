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

#define MAX_COLUMNS 32

os_mutex_t frame_mutex;
Scene      scene;
Cubemap    skybox;

int num_columns = 1;
int init_scale = 2;

bool workers_should_stop = false;
_Atomic uint32_t accum_generation = 0;
Vector3 *accum = NULL;
Vector3 *frame = NULL;
int      frame_w = 0;
int      frame_h = 0;
float accum_counts[MAX_COLUMNS] = {0};
os_mutex_t frame_mutex;
os_condvar_t accum_conds[MAX_COLUMNS];

void start_workers(os_thread *workers);
void stop_workers(os_thread *workers);

void screenshot(void);
void parse_arguments_or_exit(int argc, char **argv, int *num_columns, int *init_scale, char **scene_file);

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

Vector3 fresnelSchlick(float u, Vector3 f0)
{
	return combine(f0, combine(vec_from_scalar(1.0), f0, 1, -1), 1, pow(1.0 - u, 5.0));
}

Vector3 pixel(float x, float y, float aspect_ratio)
{
	assert(!isnan(aspect_ratio));

	Ray in_ray = ray_through_screen_at(x, y, aspect_ratio);
	assert(!isnanv(in_ray.direction));

	// Choose a light source
	int light_index = -1;
	for (int i = 0; i < scene.num_objects; i++) {
		if (scene.objects[i].material.emission_power > 0) {
			light_index = i;
			break;
		}
	}

	Vector3 contrib = {1, 1, 1};
	Vector3 result = {0, 0, 0};

	int bounces = 10;
	for (int i = 0; i < bounces; i++) {

		HitInfo hit = trace_ray(in_ray, &scene);
		if (hit.object == -1) {
			//Vector3 sky_color = {0.6, 0.7, 0.9};
			//Vector3 sky_color = {0, 0, 0};
			//Vector3 sky_color = {1, 1, 1};
			Vector3 sky_color = sample_cubemap(&skybox, normalize(in_ray.direction));
			result = combine(result, mulv(sky_color, contrib), 1, 1);
			break;
		}

		Vector3 sampled_light_color = {0, 0, 0};
		if (light_index != -1) {
			Vector3 dir_to_light_source = combine(origin_of(scene.objects[light_index]), hit.point, 1, -1);
			int max_samples = 3;
			int num_samples = 0;
			float spread = 0.5;
			for (int k = 0; k < max_samples; k++) {
				// Add some noise based on roughness
				Vector3 rand_dir = random_direction();
				if (dotv(rand_dir, hit.normal) > 0) {
					Vector3 sample_dir = normalize(combine(rand_dir, dir_to_light_source, spread, 1));
					Ray sample_ray = { combine(hit.point, sample_dir, 1, 0.001), sample_dir };
					HitInfo hit2 = trace_ray(sample_ray, &scene);
					if (hit2.object != -1)
						sampled_light_color = combine(sampled_light_color, scene.objects[hit2.object].material.emission_color, 1, scene.objects[hit2.object].material.emission_power);
					num_samples++;
				}
			}
			if (num_samples > 0)
				sampled_light_color = scale(sampled_light_color, 1.0f / num_samples);
		}

		Material material = scene.objects[hit.object].material;

		Vector3 v = scale(in_ray.direction, -1);
		Vector3 n = hit.normal;
		float NoV = clamp(dotv(n, v), 0, 1);

		Vector3 f0_dielectric = vec_from_scalar(0.16 * material.reflectance * material.reflectance);
		Vector3 f0_metal = material.albedo;
		Vector3 f0 = combine(f0_dielectric, f0_metal, (1 - material.metallic), material.metallic);
		Vector3 F = fresnelSchlick(NoV, f0);

		Vector3 rand_dir = random_direction();
		if (dotv(rand_dir, hit.normal) < 0)
			rand_dir = scale(rand_dir, -1);

		result = combine(result, mulv(scale(material.emission_color, material.emission_power), contrib), 1, 1);

		Vector3 out_dir;
		if (material.metallic > 0.001 || random_float() <= avgv(F)) {
			// Specular ray
			Vector3 reflect_dir = reflect(in_ray.direction, scale(hit.normal, -1));
			out_dir = normalize(combine(rand_dir, reflect_dir, material.roughness, 1));
		} else {
			// Diffuse ray
			out_dir = rand_dir;
			contrib = mulv(contrib, scale(material.albedo, (1 - material.metallic)));
		}
		Ray out_ray = { combine(hit.point, out_dir, 1, 0.001), out_dir };

		float light_sample_weight = 0.05;
		if (!iszerov(sampled_light_color)) {
			result = combine(result, mulv(sampled_light_color, contrib), 1, light_sample_weight);
			contrib = scale(contrib, 1 - light_sample_weight);
		}

		in_ray = out_ray;
	}

	result.x = clamp(result.x, 0, 1);
	result.y = clamp(result.y, 0, 1);
	result.z = clamp(result.z, 0, 1);

	return result;
}

float render_column(Vector3 *data, int scale_, int column_w, int column_i, int frame_w, int frame_h, uint64_t cached_generation)
{
	float scale2inv = 1.0f / (scale_ * scale_);

	int column_x = column_w * column_i;
	float aspect_ratio = (float) frame_w / frame_h;

	int lowres_frame_w = frame_w / scale_;
	int lowres_frame_h = frame_h / scale_;
	int lowres_column_w = column_w / scale_ + 1;
	int lowres_column_x = column_x / scale_;

	for (int j = 0; j < lowres_frame_h; j++) {
		for (int i = 0; i < lowres_column_w; i++) {
			float u = (float) (lowres_column_x + i) / (lowres_frame_w - 1);
			float v = (float) j  / (lowres_frame_h - 1);
			u = 1 - u;
			v = 1 - v;

			int tile_w = scale_;
			int tile_h = scale_;
			if (tile_w > column_w - i * scale_)
				tile_w = column_w - i * scale_;

			Vector3 color = pixel(u, v, aspect_ratio);
			for (int g = 0; g < tile_h; g++)
				for (int t = 0; t < tile_w; t++) {
					int pixel_index = (j * scale_ + g) * column_w + (i * scale_ + t);
					assert(pixel_index >= 0 && pixel_index < column_w * frame_h);
					data[pixel_index] = scale(color, 1);
				}
		}

		if (cached_generation != atomic_load(&accum_generation))
			break;
	}

	return scale2inv;
}

os_threadreturn worker(void *arg)
{
	float    column_data_weight = 0;
	Vector3 *column_data = NULL;
	int column_i = (int) arg;
	int column_w = 0;
	int cached_frame_w;
	int cached_frame_h;
	uint64_t cached_generation;

	int scale_ = init_scale;

	os_mutex_lock(&frame_mutex);
	while (!workers_should_stop) {
		bool resize = false;
		
		if (column_data == NULL || cached_generation != atomic_load(&accum_generation))
			resize = true;
		column_w = frame_w / num_columns;
		cached_frame_w = frame_w;
		cached_frame_h = frame_h;
		cached_generation = atomic_load(&accum_generation);
		os_mutex_unlock(&frame_mutex);

		if (resize) {
			free(column_data);
			column_data = malloc(sizeof(Vector3) * column_w * cached_frame_h);
			if (!column_data) abort();
		}

		column_data_weight += render_column(column_data, scale_, column_w, column_i, cached_frame_w, cached_frame_h, cached_generation);

		os_mutex_lock(&frame_mutex);

		// Publish changes
		if (cached_generation == atomic_load(&accum_generation)) {
			for (int j = 0; j < frame_h; j++)
				for (int i = 0; i < column_w; i++) {	
					int column_x = column_w * column_i;
					int src_index = j * column_w + i;
					int dst_index = j * frame_w + (i + column_x);
					assert(src_index >= 0 && src_index < column_w * cached_frame_h);
					assert(dst_index >= 0 && dst_index < cached_frame_w * cached_frame_h);
					accum[dst_index] = combine(accum[dst_index], column_data[src_index], 1, 1.0f / (scale_ * scale_));
				}
			os_condvar_signal(&accum_conds[column_i]);
			accum_counts[column_i] += column_data_weight;
			if (scale_ > 1)
				scale_ >>= 1;
		} else {
			scale_ = init_scale;
		}

		column_data_weight = 0;
	}
	os_mutex_unlock(&frame_mutex);
}

void update_frame_texture(void)
{
	os_mutex_lock(&frame_mutex);

	if (frame_w != get_screen_w() || frame_h != get_screen_h()) {

		frame_w = get_screen_w();
		frame_h = get_screen_h();

		if (frame) free(frame);
		if (accum) free(accum);

		frame = malloc(sizeof(Vector3) * frame_w * frame_h);
		if (!frame) { printf("OUT OF MEMORY\n"); abort(); }

		accum = malloc(sizeof(Vector3) * frame_w * frame_h);
		if (!accum) { printf("OUT OF MEMORY\n"); abort(); }

		for (int i = 0; i < num_columns; i++)
			accum_counts[i] = 0;
		
		memset(accum, 0, sizeof(Vector3) * frame_w * frame_h);
		memset(frame, 0, sizeof(Vector3) * frame_w * frame_h);

		atomic_fetch_add(&accum_generation, 1);
	}

	int column_w = frame_w / num_columns;

	for (int i = 0; i < num_columns; i++) {
		while (accum_counts[i] < 0.0001)
			os_condvar_wait(&accum_conds[i], &frame_mutex, -1);
	}

	for (int j = 0; j < frame_h; j++)
		for (int i = 0; i < frame_w; i++) {

			float u = (float) i / (frame_w - 1);
			float v = (float) j / (frame_h - 1);
			u = 1 - u;
			v = 1 - v;

			int pixel_index = j * frame_w + i;
			frame[pixel_index] = scale(accum[pixel_index], 1.0f / accum_counts[i / column_w]);
		}

	move_frame_to_the_gpu(frame_w, frame_h, frame);
	os_mutex_unlock(&frame_mutex);
}

int main(int argc, char **argv)
{
	char *scene_file;
	parse_arguments_or_exit(argc, argv, &num_columns, &init_scale, &scene_file);

	if (!parse_scene_file(scene_file, &scene))
		return -1;

	const char *faces[] = {
		[CF_RIGHT]  = "assets/skybox/right.jpg",
		[CF_LEFT]   = "assets/skybox/left.jpg",
		[CF_TOP]    = "assets/skybox/top.jpg",
		[CF_BOTTOM] = "assets/skybox/bottom.jpg",
		[CF_FRONT]  = "assets/skybox/front.jpg",
		[CF_BACK]   = "assets/skybox/back.jpg",
	};
	load_cubemap(&skybox, faces);

	startup_window_and_opengl_context_or_exit(2 * 640, 2 * 480, "Ray Tracing");

	os_thread workers[MAX_COLUMNS];
	start_workers(workers);

	for (;;) {

		bool exit = false;
		for (;;) {

			double mouse_x;
			double mouse_y;
			int event = pop_event(&mouse_x, &mouse_y);

			float speed = 0.5;
			if (event == EVENT_CLOSE || event == EVENT_PRESS_ESC) {

				exit = true;
				break;

			} else if (event == EVENT_PRESS_W || event == EVENT_AGAIN_W) {

				move_camera(UP, speed);
				invalidate_accumulation();

			} else if (event == EVENT_PRESS_A || event == EVENT_AGAIN_A) {

				move_camera(LEFT, speed);
				invalidate_accumulation();

			} else if (event == EVENT_PRESS_S || event == EVENT_AGAIN_S) {

				move_camera(DOWN, speed);
				invalidate_accumulation();

			} else if (event == EVENT_PRESS_D || event == EVENT_AGAIN_D) {

				move_camera(RIGHT, speed);
				invalidate_accumulation();

			} else if (event == EVENT_MOVE_MOUSE) {

				rotate_camera(mouse_x, mouse_y);
				invalidate_accumulation();

			} else if (event == EVENT_PRESS_SPACE) {
				screenshot();
			}
		}
		if (exit) break;

		update_frame_texture();
		draw_frame();
	}

	stop_workers(workers);
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

// Must be executed on the main thread
void screenshot(void)
{
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

void start_workers(os_thread *workers)
{
	os_mutex_create(&frame_mutex);

	for (int i = 0; i < num_columns; i++)
		os_condvar_create(&accum_conds[i]);

	for (int i = 0; i < num_columns; i++)
		os_thread_create(&workers[i], (void*) i, worker);
}

void stop_workers(os_thread *workers)
{
	os_mutex_lock(&frame_mutex);
	workers_should_stop = true;
	os_mutex_unlock(&frame_mutex);
	for (int i = 0; i < num_columns; i++)
		os_thread_join(workers[i]);

	for (int i = 0; i < num_columns; i++)
		os_condvar_delete(&accum_conds[i]);
}
