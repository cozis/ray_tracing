#include <stdio.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include "utils.h"
#include "gpu_and_windowing.h"

static GLFWwindow *window;
static unsigned int screen_program;
static unsigned int frame_texture;
static unsigned int vao;
static unsigned int vbo;
static int screen_w;
static int screen_h;

#define MAX_EVENTS 512
int event_queue[MAX_EVENTS];
int event_queue_head = 0;
int event_queue_size = 0;

void load_cubemap(Cubemap *c, const char *files[6])
{
	for (int i = 0; i < 6; i++) {
		c->data[i] = stbi_load(files[i], &c->w, &c->h, &c->chan, 0);
		if (c->data[i] == NULL) {
			fprintf(stderr, "Couldn't load image '%s'\n", files[i]);
			abort();
		}
	}
}

void free_cubemap(Cubemap *c)
{
	for (int i = 0; i < 6; i++) {
		stbi_image_free(c->data[i]);
	}
}

Vector3 sample_cubemap(Cubemap *c, Vector3 dir)
{
	float abs_x = absf(dir.x);
	float abs_y = absf(dir.y);
	float abs_z = absf(dir.z);

	CubeFace face;

	float u;
	float v;
	float eps = 0;

	if (abs_x > abs_y && abs_x > abs_z) {
		// X dominant
		if (dir.x > 0) {
			// right face
			face = CF_RIGHT;
			u = -dir.z / (abs_x + eps);
			v = -dir.y / (abs_x + eps);
		} else {
			// left face
			face = CF_LEFT;
			u = dir.z / (abs_x + eps);
			v = -dir.y / (abs_x + eps);
		}
	} else if (abs_y > abs_x && abs_y > abs_z) {
		// Y dominant
		assert(abs_y > 0);
		if (dir.y > 0) {
			// top face
			face = CF_TOP;
			u = dir.x / (abs_y + eps);
			v = dir.z / (abs_y + eps);
		} else {
			// bottom face
			face = CF_BOTTOM;
			u = dir.x / (abs_y + eps);
			v = -dir.z / (abs_y + eps);
		}
	} else {
		// Z dominant
		if (dir.z > 0) {
			// front face
			face = CF_FRONT;
			u = dir.x / (abs_z + eps);
			v = -dir.y / (abs_z + eps);
		} else {
			// back face
			face = CF_BACK;
			u = -dir.x / (abs_z + eps);
			v = -dir.y / (abs_z + eps);
		}
	}

	u = clamp(u, -1, 1);
	v = clamp(v, -1, 1);

	u = 0.5f * (u + 1.0f);
	v = 0.5f * (v + 1.0f);

	// Pixel coordinates
	int x = u * (c->w - 1);
	int y = v * (c->h - 1);

	uint8_t *color = &c->data[face][(y * c->w + x) * c->chan];
	return (Vector3) {
		(float) color[0] / 255,
		(float) color[1] / 255,
		(float) color[2] / 255,
	};
}

static unsigned int
compile_shader(const char *vertex_file, const char *fragment_file)
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
	glShaderSource(vertex_shader, 1, (const GLchar * const *) &vertex_str, NULL);
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
	glShaderSource(fragment_shader, 1, (const GLchar * const *) &fragment_str, NULL);
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

static void push_event(int event)
{
	if (event_queue_size == MAX_EVENTS) {
		fprintf(stderr, "Event queue full. An event has been lost\n");
		return;
	}
	int tail = (event_queue_head + event_queue_size) % MAX_EVENTS;
	event_queue[tail] = event;
	event_queue_size++;
}

int pop_event(double *mouse_x, double *mouse_y)
{
	if (glfwWindowShouldClose(window))
		return EVENT_CLOSE;

	if (event_queue_size == 0)
		return EVENT_EMPTY;

	int event = event_queue[event_queue_head];
	event_queue_head = (event_queue_head + 1) % MAX_EVENTS;
	event_queue_size--;

	if (event == EVENT_MOVE_MOUSE)
		glfwGetCursorPos(window, mouse_x, mouse_y);
	return event;
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	switch (key) {
		case GLFW_KEY_SPACE:
		if (action == GLFW_PRESS) push_event(EVENT_PRESS_SPACE);
		break;

		case GLFW_KEY_ESCAPE:
		if (action == GLFW_PRESS) push_event(EVENT_PRESS_ESC);
		break;
	}
}

static void cursor_callback(GLFWwindow *window, double x, double y)
{
	push_event(EVENT_MOVE_MOUSE);
}

static void framebuffer_size_callback(GLFWwindow* window, int w, int h)
{
    glViewport(0, 0, w, h);
	screen_w = w;
	screen_h = h;
}

void startup_window_and_opengl_context_or_exit(int window_w, int window_h, const char *title)
{
	glfwSetErrorCallback(error_callback);

	if (!glfwInit())
		exit(-1);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(window_w, window_h, title, NULL, NULL);
	if (!window) {
		glfwTerminate();
		exit(-1);
	}

	glfwSetKeyCallback(window, key_callback);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		printf("Failed to initialize GLAD\n");
		exit(-1);
	}

	glfwSwapInterval(1);

	glfwGetWindowSize(window, &screen_w, &screen_h);

	screen_program = compile_shader("assets/screen.vs", "assets/screen.fs");
	if (!screen_program) {
		printf("Couldn't compile program\n");
		exit(-1);
	}
	set_uniform_i(screen_program, "screenTexture", 0);

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
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void cleanup_window_and_opengl_context(void)
{
	glfwDestroyWindow(window);
	glfwTerminate();
}

int get_screen_w(void)
{
	return screen_w;
}

int get_screen_h(void)
{
	return screen_h;
}

void move_frame_to_the_gpu(int w, int h, Vector3 *data)
{
	glBindTexture(GL_TEXTURE_2D, frame_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_FLOAT, data);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void draw_frame(void)
{
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	glUseProgram(screen_program);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, frame_texture);
	glBindVertexArray(vao);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);

	glfwSwapBuffers(window);
	glfwPollEvents();

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) { push_event(EVENT_PRESS_W); }
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) { push_event(EVENT_PRESS_A); }
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) { push_event(EVENT_PRESS_S); }
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { push_event(EVENT_PRESS_D); }
}