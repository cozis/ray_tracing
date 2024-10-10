#include <stdint.h>
#include "vector.h"

typedef struct {
	uint8_t *data[6];
	int w, h, chan;
} Cubemap;

typedef enum {
	CF_FRONT,
	CF_BACK,
	CF_LEFT,
	CF_RIGHT,
	CF_TOP,
	CF_BOTTOM,
} CubeFace;

enum {
	EVENT_EMPTY = 0,
	EVENT_CLOSE,
	EVENT_PRESS_SPACE,
	EVENT_PRESS_ESC,
	EVENT_PRESS_W,
	EVENT_PRESS_A,
	EVENT_PRESS_S,
	EVENT_PRESS_D,
	EVENT_AGAIN_SPACE,
	EVENT_AGAIN_ESC,
	EVENT_AGAIN_W,
	EVENT_AGAIN_A,
	EVENT_AGAIN_S,
	EVENT_AGAIN_D,
	EVENT_MOVE_MOUSE,
};

int pop_event(double *mouse_x, double *mouse_y);

void startup_window_and_opengl_context_or_exit(int window_w, int window_h, const char *title);
void cleanup_window_and_opengl_context(void);

int get_screen_w(void);
int get_screen_h(void);

void move_frame_to_the_gpu(int w, int h, Vector3 *data);
void draw_frame(void);

void    load_cubemap(Cubemap *c, const char *files[6]);
void    free_cubemap(Cubemap *c);
Vector3 sample_cubemap(Cubemap *c, Vector3 dir);
