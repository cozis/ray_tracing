#include <math.h>
#include "camera.h"

static bool first_mouse = true;
static float yaw = -90.0f;
static float pitch = 0.0f;
static float last_x = 800.0f / 2.0;
static float last_y = 600.0f / 2.0;
static float fov    = 30.0f;

static float delta_time = 0.0f;
static float last_frame = 0.0f;

static Vector3 camera_pos   = {5, 5, 5};
static Vector3 camera_front = {-1, -1, -1};
static Vector3 camera_up    = {0, 1, 0};

Vector3 get_camera_pos(void)
{
	return camera_pos;
}

void rotate_camera(double mouse_x, double mouse_y)
{
	float x = mouse_x;
	float y = mouse_y;

	if (first_mouse) {
		last_x = x;
		last_y = y;
		first_mouse = false;
	}

	float dx = x - last_x;
	float dy = last_y - y;
	last_x = x;
	last_y = y;

	float sensitivity = 0.1f;
	dx *= sensitivity;
	dy *= sensitivity;

	yaw   += dx;
	pitch += dy;

	if (pitch >  89.0f) pitch =  89.0f;
	if (pitch < -89.0f) pitch = -89.0f;

	float   yaw_rad = deg2rad(yaw);
	float pitch_rad = deg2rad(pitch);

	Vector3 front;
	front.x = cos(yaw_rad) * cos(pitch_rad);
	front.y = sin(pitch_rad);
	front.z = sin(yaw_rad) * cos(pitch_rad);
	front = normalize(front);
	
	camera_front = front;
}

void move_camera(Direction dir, float speed)
{
    switch (dir) {
        case UP   : camera_pos = combine(camera_pos, camera_front, 1, +speed); break;
        case DOWN : camera_pos = combine(camera_pos, camera_front, 1, -speed); break;
        case LEFT : camera_pos = combine(camera_pos, normalize(cross(camera_front, camera_up)), 1, -speed); break;
        case RIGHT: camera_pos = combine(camera_pos, normalize(cross(camera_front, camera_up)), 1, +speed); break;
    }
}

Matrix4 camera_pov(void)
{
	return lookat_matrix(camera_pos, combine(camera_pos, camera_front, 1, 1), camera_up);
}

Ray ray_through_screen_at(float px, float py, float aspect_ratio)
{
	Vector3 w = normalize(scale(camera_front, -1));
	Vector3 u = normalize(cross(camera_up, w));
	Vector3 v = cross(w, u);

	float screen_h = 2 * tan(fov / 2);
	float screen_w = aspect_ratio * screen_h;

	Vector3 horizontal = scale(u, screen_w);
	Vector3 vertical   = scale(v, screen_h);

	Vector3 lower_left_corner = combine4(camera_pos, horizontal, vertical, w, 1, -0.5, -0.5, -1);

	Vector3 dir = combine4(lower_left_corner, horizontal, vertical, camera_pos, 1, px, py, -1);

	return (Ray) {.origin=camera_pos, .direction=dir};
}
