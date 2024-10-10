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
#include <assert.h>
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
	assert(!isnan(aspect_ratio));

	Vector3 w = normalize(scalev(camera_front, -1));
	Vector3 u = normalize(cross(camera_up, w));
	Vector3 v = cross(w, u);

	assert(!isnanv(w));
	assert(!isnanv(u));
	assert(!isnanv(v));

	float screen_h = 2 * tan(fov / 2);
	float screen_w = aspect_ratio * screen_h;
	assert(!isnan(screen_h));
	assert(!isnan(screen_w));

	Vector3 horizontal = scalev(u, screen_w);
	Vector3 vertical   = scalev(v, screen_h);

	assert(!isnanv(horizontal));
	assert(!isnanv(vertical));

	Vector3 lower_left_corner = combine4(camera_pos, horizontal, vertical, w, 1, -0.5, -0.5, -1);
	assert(!isnanv(lower_left_corner));

	Vector3 dir = combine4(lower_left_corner, horizontal, vertical, camera_pos, 1, px, py, -1);
	assert(!isnanv(dir));

	return (Ray) {.origin=camera_pos, .direction=dir};
}
