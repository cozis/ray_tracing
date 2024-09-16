#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h> // FLT_MAX
#include <glad/glad.h>
//#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "utils.h"
#include "camera.h"
#include "vector.h"
#include "mesh.h"

typedef struct {
	Vector3 albedo;
	float   metallic;
	float   roughness;
	float   emission_power;
	Vector3 emission_color;
} Material;

int screen_w;
int screen_h;

float maxf(float x, float y) { return x > y ? x : y; }
float minf(float x, float y) { return x < y ? x : y; }

static unsigned int
compile_shader(const char *vertex_file,
               const char *fragment_file)
{
    int  success;
    char infolog[512];

    char *vertex_str = load_file(vertex_file, NULL);
    if (vertex_str == NULL) {
        fprintf(stderr, "Couldn't load file '%s'\n", vertex_file);
        return 0;
    }

    char *fragment_str = load_file(fragment_file, NULL);
    if (fragment_str == NULL) {
        fprintf(stderr, "Couldn't load file '%s'\n", fragment_file);
        free(vertex_str);
        return 0;
    }

    unsigned int vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_str, NULL);
    glCompileShader(vertex_shader);

    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if(!success) {
        glGetShaderInfoLog(vertex_shader, sizeof(infolog), NULL, infolog);
        fprintf(stderr, "Couldn't compile vertex shader '%s' (%s)\n", vertex_file, infolog);
        free(vertex_str);
        free(fragment_str);
        return 0;
    }

    unsigned int fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_str, NULL);
    glCompileShader(fragment_shader);

    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if(!success) {
        glGetShaderInfoLog(fragment_shader, sizeof(infolog), NULL, infolog);
        fprintf(stderr, "Couldn't compile fragment shader '%s' (%s)\n", fragment_file, infolog);
        free(vertex_str);
        free(fragment_str);
        return 0;
    }

    unsigned int shader_program = glCreateProgram();
    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    glLinkProgram(shader_program);

    glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(shader_program, sizeof(infolog), NULL, infolog);
        fprintf(stderr, "Couldn't link shader program (%s)\n", infolog);
        free(vertex_str);
        free(fragment_str);
        return 0;
    }

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    free(vertex_str);
    free(fragment_str);
    return shader_program;
}

static void set_uniform_m4(unsigned int program, const char *name, Matrix4 value)
{
	int location = glGetUniformLocation(program, name);
	if (location < 0) {
		printf("Can't set uniform '%s'\n", name);
		abort();
	}
	glUniformMatrix4fv(location, 1, GL_FALSE, (float*) &value);
}

static void set_uniform_v3(unsigned int program, const char *name, Vector3 value)
{
	int location = glGetUniformLocation(program, name);
	if (location < 0) {
		printf("Can't set uniform '%s' (program %d, location %d)\n", name, program, location);
		abort();
	}
	glUniform3f(location, value.x, value.y, value.z);
}

static void set_uniform_i(unsigned int program, const char *name, int value)
{
	int location = glGetUniformLocation(program, name);
	if (location < 0) {
		printf("Can't set uniform '%s'\n", name);
		abort();
	}
	glUniform1i(location, value);
}

static void set_uniform_f(unsigned int program, const char *name, float value)
{
	int location = glGetUniformLocation(program, name);
	if (location < 0) {
		printf("Can't set uniform '%s'\n", name);
		abort();
	}
	glUniform1f(location, value);
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void reset_accum(void);

void cursor_pos_callback(GLFWwindow *window, double x, double y)
{
	reset_accum();
    rotate_camera(x, y);
}

typedef struct {
	Vector3 origin;
	Vector3 size;
} Cube;

bool intersect_cube(Ray r, Cube c, float *tnear, float *tfar, Vector3 *normal)
{
	float txmin, txmax;
	float tymin, tymax;
	float tzmin, tzmax;

	float tn;
	float tf;

	Vector3 a = c.origin;
	Vector3 b = combine(c.origin, c.size, 1, 1);

	int hit_axis = 0; // 0=x, 1=y, 2=z

	if (r.direction.x >= 0) {
		txmin = (a.x - r.origin.x) / r.direction.x;
		txmax = (b.x - r.origin.x) / r.direction.x;
	} else {
		txmax = (a.x - r.origin.x) / r.direction.x;
		txmin = (b.x - r.origin.x) / r.direction.x;
	}

	if (r.direction.y >= 0) {
		tymin = (a.y - r.origin.y) / r.direction.y;
		tymax = (b.y - r.origin.y) / r.direction.y;
	} else {
		tymax = (a.y - r.origin.y) / r.direction.y;
		tymin = (b.y - r.origin.y) / r.direction.y;
	}

	if (txmin > tymax || tymin > txmax)
		return false;

	if (tymin > txmin) { txmin = tymin; hit_axis = 1; }
	if (tymax < txmax) txmax = tymax;

	if (r.direction.z >= 0) {
		tzmin = (a.z - r.origin.z) / r.direction.z;
		tzmax = (b.z - r.origin.z) / r.direction.z;
	} else {
		tzmax = (a.z - r.origin.z) / r.direction.z;
		tzmin = (b.z - r.origin.z) / r.direction.z;
	}

	if (txmin > tzmax || tzmin > txmax)
		return false;
	
	if (tzmin > txmin) { txmin = tzmin; hit_axis = 2; };
	if (tzmax < txmax) txmax = tzmax;

	if (tnear) *tnear = txmin;
	if (tfar)  *tfar  = txmax;
	if (normal) {
		switch (hit_axis) {
			case 0: *normal = r.direction.x > 0 ? (Vector3) {-1, 0, 0} : (Vector3) {1, 0, 0}; break;
			case 1: *normal = r.direction.y > 0 ? (Vector3) {0, -1, 0} : (Vector3) {0, 1, 0}; break;
			case 2: *normal = r.direction.z > 0 ? (Vector3) {0, 0, -1} : (Vector3) {0, 0, 1}; break;
		}
	}
	return true;
}

bool intersect_sphere(Ray r, Sphere s, float *t)
{
	/*
	 * Any point of the ray can be written as
	 *
	 *     P(t) = O + t * D
	 * 
	 * with O origin and D direction.
	 * 
	 * All points P=(x,y,z) of a sphere can be described as
	 * those (and only those) that satisfy the equation
	 *
	 *     x^2 + y^2 + z^2 = R^2
	 *     P^2 - R^2 = 0
	 * 
	 * with R radius of the sphere. The sphere here is centered
	 * at the origin.
	 * 
	 * Intersection points of the ray with the sphere must satisfy
	 * both:
	 * 
	 *     P(t) = O + t * D
	 *     P^2 - R^2 = 0
	 * 
	 *     => (O + tD)^2 - R^2 = 0
	 *     => t^2 * D^2 + t * 2OD + O^2 - R^2 = 0
	 * 
	 * we can use the quadratic formula here, and more specifically
	 * the discriminant to check if solutions exist and how many
	 */
	Vector3 oc = combine(s.center, r.origin, 1, -1);
	float a = dotv(r.direction, r.direction);
	float b = -2 * dotv(oc, r.direction);
	float c = dotv(oc, oc) - s.radius * s.radius;

	float discr = b*b - 4*a*c;

	if (discr > 0) {
		float s0 = (- b + sqrt(discr)) / (2 * a);
		float s1 = (- b - sqrt(discr)) / (2 * a);
		if (s0 > s1) {
			float tmp = s0;
			s0 = s1;
			s1 = tmp;
		}
		if (s0 < 0) {
			s0 = s1;
			if (s0 < 0) return false;
		}
		if (t) *t = s0;
		return true;
	}

	// Zero solutions
	return false;
}

typedef enum {
	OBJECT_CUBE,
	OBJECT_SPHERE,
} ObjectType;

typedef struct {
	ObjectType type;
	union {
		Sphere sphere;
		Cube cube;
	};
	Material material;
} Object;

Object   cube(Material material, Vector3 origin, Vector3 size) { return (Object) {.material=material, .type=OBJECT_CUBE, .cube=(Cube) {.origin=origin, .size=size}}; }
Object sphere(Material material, Vector3 origin, float radius) { return (Object) {.material=material, .type=OBJECT_SPHERE, .sphere=(Sphere) {.center=origin, .radius=radius}}; }

bool intersect_object(Ray r, Object o, float *t, Vector3 *normal)
{
	switch (o.type) {

		case OBJECT_CUBE:
		return intersect_cube(r, o.cube, t, NULL, normal);

		case OBJECT_SPHERE:
		if (intersect_sphere(r, o.sphere, t)) {
			if (normal) {
				Vector3 hit_point = combine(r.origin, r.direction, 1, *t);
				*normal = normalize(combine(hit_point, o.sphere.center, 1, -1));
			}
			return true;
		}
		return false;
	}
	return false;
}

unsigned int wang_hash(unsigned int seed)
{
	seed = (seed ^ 61) ^ (seed >> 16);
	seed *= 9;
	seed = seed ^ (seed >> 4);
	seed *= 0x27d4eb2d;
	seed = seed ^ (seed >> 15);
	return seed;
}

unsigned int pcg_hash(unsigned int input)
{
	unsigned int state = input * 747796405U + 2891336453U;
	unsigned int word = ((state >> ((state >> 28U) + 4U)) ^ state) * 277803737U;
	return (word >> 22U) ^ word;
}

float random_float(unsigned int *seed)
{
//	*seed = pcg_hash(*seed);
//	return (float) *seed / UINT_MAX;
	return (float) rand() / RAND_MAX;
}

Vector3 random_vector(unsigned int seed)
{
	return (Vector3) {
		.x = random_float(&seed) * 2 - 1,
		.y = random_float(&seed) * 2 - 1,
		.z = random_float(&seed) * 2 - 1,
	};
}

Vector3 random_direction(unsigned int seed)
{
	return normalize(random_vector(seed));
}

Vector3 reflect(Vector3 dir, Vector3 normal)
{
	float f = -2 * dotv(normal, dir);
	return combine(dir, normal, 1, f);
}

#define MAX_OBJECTS 32
Object objects[MAX_OBJECTS];
int num_objects = 0;

void add_object(Object o)
{
	if (num_objects < MAX_OBJECTS)
		objects[num_objects++] = o;
}

typedef struct {
	float   distance;
	Vector3 point;
	Vector3 normal;
	int     object;
} HitInfo;

HitInfo trace_ray(Ray ray)
{
	ray.direction = normalize(ray.direction);

	float   nearest_t = FLT_MAX;
	int     nearest_object = -1;
	Vector3 nearest_normal;
	for (int i = 0; i < num_objects; i++) {
		float t;
		Vector3 n;
		if (!intersect_object(ray, objects[i], &t, &n))
			continue;
		if (t >= 0 && t < nearest_t) {
			nearest_t = t;
			nearest_object = i;
			nearest_normal = n;
		}
	}

	if (nearest_object == -1) {
		HitInfo result;
		result.distance = -1;
		result.normal   = (Vector3) {0, 0, 0};
		result.point    = (Vector3) {0, 0, 0};
		result.object   = -1;
		return result;
	} else {
		HitInfo result;
		result.distance = nearest_t;
		result.normal   = nearest_normal;
		result.point    = combine(ray.origin, ray.direction, 1, nearest_t);
		result.object   = nearest_object;
		return result;
	}
}

Vector3 pixel(float x, float y)
{
	Ray ray = ray_through_screen_at(x, y, (float) screen_w/screen_h);

	Vector3 contrib = {1, 1, 1};
	Vector3 light = {0, 0, 0};

	int rays_per_pixel = 1;
	for (int j = 0; j < rays_per_pixel; j++) {

		int bounces = 4;
		for (int i = 0; i < bounces; i++) {

			HitInfo hit = trace_ray(ray);
			if (hit.object == -1) {
				//Vector3 sky_color = {0.6, 0.7, 0.9};
				Vector3 sky_color = {0, 0, 0};
				light = combine(light, mulv(sky_color, contrib), 1, 1);
				break;
			}

			Material material = objects[hit.object].material;
			contrib = mulv(contrib, material.albedo);
			light = combine(light, material.emission_color, 1, material.emission_power);
#if 0
			Vector3 reflect_dir = reflect(ray.direction, scale(hit.normal, -1));
			Vector3 noise_dir = scale(random_direction(), 0.5);
			if (dotv(noise_dir, reflect_dir) < 0)
				noise_dir = scale(noise_dir, -1);

			float roughness = objects[hit.object].material.roughness;
			Vector3 new_dir = combine(noise_dir, reflect_dir, roughness, 1);
#endif

			Vector3 new_dir = random_direction(i * 1000000 + x * 1000 + y);
/*
			if (dotv(new_dir, hit.normal) < 0)
				new_dir = scale(new_dir, -1);
*/
			ray = (Ray) { combine(hit.point, new_dir, 1, 0.001), new_dir };
		}
	}
	light = scale(light, 1.0f/rays_per_pixel);

	return light;
}

Vector3 *accum = NULL;
Vector3 *frame = NULL;
int      frame_w = 0;
int      frame_h = 0;
unsigned int frame_texture;
int      accum_index = 1;

void reset_accum(void)
{
	accum_index = 1;
}

void update_frame_texture(float s)
{
	if (frame_w != s * screen_w || frame_h != s * screen_h) {
		frame_w = s * screen_w;
		frame_h = s * screen_h;

		if (frame) free(frame);
		if (accum) free(accum);

		frame = malloc(sizeof(Vector3) * frame_w * frame_h);
		if (!frame) { printf("OUT OF MEMORY\n"); abort(); }

		accum = malloc(sizeof(Vector3) * frame_w * frame_h);
		if (!accum) { printf("OUT OF MEMORY\n"); abort(); }

		accum_index = 1;
	}

	if (accum_index == 1)
		memset(accum, 0, sizeof(Vector3) * frame_w * frame_h);

	for (int j = 0; j < frame_h; j++)
		for (int i = 0; i < frame_w; i++) {
			float u = (float) i / (frame_w - 1);
			float v = (float) j / (frame_h - 1);
			u = 1 - u;
			v = 1 - v;

			Vector3 color = pixel(u, v);

			int pixel_index = j * frame_w + i;
			accum[pixel_index] = combine(accum[pixel_index], color, 1, 1);
			frame[pixel_index] = scale(accum[pixel_index], 1.0f / accum_index);
		}
	
	accum_index++;

	glBindTexture(GL_TEXTURE_2D, frame_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frame_w, frame_h, 0, GL_RGB, GL_FLOAT, frame);
	glBindTexture(GL_TEXTURE_2D, 0);
}

int main(void)
{
/*
	add_object(cube((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=1, .albedo=(Vector3) {1, 0, 0}}, (Vector3) {0, 0, 0}, (Vector3) {10, 5, 0.1})),
	add_object(cube((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=1, .albedo=(Vector3) {1, 0, 0}}, (Vector3) {0, 0, 0}, (Vector3) {0.1, 5, 10})),
*/
	add_object(cube((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=1, .albedo=(Vector3) {0.4, 0.3, 0.9}}, (Vector3) {0, -0.1, 0}, (Vector3) {10, 0.1, 10})),
	add_object(cube((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=1, .albedo=(Vector3) {1, 0, 0}},     (Vector3) {7, 0, 8}, (Vector3) {1, 1, 1})),
	add_object(cube((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=1, .albedo=(Vector3) {1, 0, 1}},   (Vector3) {6, 0, 7}, (Vector3) {1, 1, 1})),
	add_object(sphere((Material) {.emission_color={1, 0, 0}, .emission_power=0.3, .metallic=0, .roughness=0, .albedo=(Vector3) {1, 0, 0}},   (Vector3) {3, 1, 3}, 1)),
	add_object(sphere((Material) {.emission_color={0}, .emission_power=0, .metallic=0, .roughness=0, .albedo=(Vector3) {0, 1, 0}},   (Vector3) {5, 1, 3}, 1)),

    glfwSetErrorCallback(error_callback);

    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow *window = glfwCreateWindow(640, 480, "Path Trace", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

	glfwGetWindowSize(window, &screen_w, &screen_h);

    glfwSetKeyCallback(window, key_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        printf("Failed to initialize GLAD\n");
        return -1;
    }

    glfwSwapInterval(1);

	unsigned int screen_program = compile_shader("assets/screen.vs", "assets/screen.fs");
	if (!screen_program) { printf("Couldn't compile program\n"); return -1; }
	set_uniform_i(screen_program, "screenTexture", 0);

	unsigned int vao, vbo;
	{
		float vertices[] = {
			// positions   // texCoords
			-1.0f,  1.0f,  0.0f, 1.0f,
			-1.0f, -1.0f,  0.0f, 0.0f,
			1.0f, -1.0f,   1.0f, 0.0f,

			-1.0f,  1.0f,  0.0f, 1.0f,
			1.0f, -1.0f,   1.0f, 0.0f,
			1.0f,  1.0f,   1.0f, 1.0f
		};

		glGenVertexArrays(1, &vao);
		glGenBuffers(1, &vbo);

		glBindVertexArray(vao);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &vertices, GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);

		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
	}

	glGenTextures(1, &frame_texture);
	glBindTexture(GL_TEXTURE_2D, frame_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	while (!glfwWindowShouldClose(window)) {

		glfwGetWindowSize(window, &screen_w, &screen_h);

		float speed = 0.5;
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) { move_camera(UP, speed); accum_index = 1; }
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) { move_camera(DOWN, speed); accum_index = 1; }
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) { move_camera(LEFT, speed); accum_index = 1; }
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { move_camera(RIGHT, speed); accum_index = 1; }

		Vector3 clear_color = {1, 1, 1};

		update_frame_texture(1);

		glViewport(0, 0, screen_w, screen_h);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, 1.0f);
		glClearStencil(0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

		glUseProgram(screen_program);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, frame_texture);
		glBindVertexArray(vao);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
