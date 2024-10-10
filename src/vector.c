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
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "vector.h"

#define EPSILON 0.00001

float maxf(float x, float y)
{
	return x > y ? x : y;
}

float minf(float x, float y)
{
	return x < y ? x : y;
}

float absf(float x)
{
	return x < 0 ? -x : x;
}

float clamp(float x, float min, float max)
{
	assert(min <= max);
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

Vector3 maxv(Vector3 a, Vector3 b)
{
	return (Vector3) {
		maxf(a.x, b.x),
		maxf(a.y, b.y),
		maxf(a.z, b.z),
	};
}

Vector3 vec_from_scalar(float s)
{
	return (Vector3) {s, s, s};
}

bool isnanv(Vector3 v)
{
	return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

bool iszerof(float f)
{
	return f < 0.0001 && f > -0.0001;
}

bool iszerov(Vector3 v)
{
	return iszerof(v.x) && iszerof(v.y) && iszerof(v.z);
}

float avgv(Vector3 v)
{
	return (v.x + v.y + v.z) / 3;
}

float deg2rad(float deg)
{
	return 3.14159265358979323846 * deg / 180;
}

Vector3 random_vector(void)
{
	return (Vector3) {
		.x = random_float() * 2 - 1,
		.y = random_float() * 2 - 1,
		.z = random_float() * 2 - 1,
	};
}

Vector3 random_direction(void)
{
	return normalize(random_vector());
}

Vector3 reflect(Vector3 dir, Vector3 normal)
{
	float f = -2 * dotv(normal, dir);
	return combine(dir, normal, 1, f);
}

float norm2_of(Vector3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

float norm_of(Vector3 v)
{
	return sqrt(norm2_of(v));
}

Vector3 normalize(Vector3 v)
{
	float norm = norm_of(v);
	if (norm < EPSILON && norm > -EPSILON)
		return v;
	v.x /= norm;
	v.y /= norm;
	v.z /= norm;
	return v;
}

Vector3 scale(Vector3 v, float f)
{
	v.x *= f;
	v.y *= f;
	v.z *= f;
	return v;
}

Vector3 combine(Vector3 u, Vector3 v, float a, float b)
{
	Vector3 r;
	r.x = u.x * a + v.x * b;
	r.y = u.y * a + v.y * b;
	r.z = u.z * a + v.z * b;
	return r;
}

Vector3 combine4(Vector3 u, Vector3 v, Vector3 g, Vector3 t, float a, float b, float c, float d)
{
	Vector3 r;
	r.x = u.x * a + v.x * b + g.x * c + t.x * d;
	r.y = u.y * a + v.y * b + g.y * c + t.y * d;
	r.z = u.z * a + v.z * b + g.z * c + t.z * d;
	return r;
}

Vector3 cross(Vector3 u, Vector3 v)
{
	Vector3 r;
	r.x = u.y * v.z - u.z * v.y;
	r.y = u.z * v.x - u.x * v.z;
	r.z = u.x * v.y - u.y * v.x;
	return r;
}

Matrix4 translate_matrix(Vector3 v, float f)
{
	Matrix4 m;
	m.data[0][0] = 1;
	m.data[0][1] = 0;
	m.data[0][2] = 0;
	m.data[0][3] = 0;
	m.data[1][0] = 0;
	m.data[1][1] = 1;
	m.data[1][2] = 0;
	m.data[1][3] = 0;
	m.data[2][0] = 0;
	m.data[2][1] = 0;
	m.data[2][2] = 1;
	m.data[2][3] = 0;
	m.data[3][0] = f * v.x;
	m.data[3][1] = f * v.y;
	m.data[3][2] = f * v.z;
	m.data[3][3] = 1;
	return m;
}

Matrix4 scale_matrix(Vector3 v)
{
	Matrix4 m;
	m.data[0][0] = v.x;
	m.data[0][1] = 0;
	m.data[0][2] = 0;
	m.data[0][3] = 0;
	m.data[1][0] = 0;
	m.data[1][1] = v.y;
	m.data[1][2] = 0;
	m.data[1][3] = 0;
	m.data[2][0] = 0;
	m.data[2][1] = 0;
	m.data[2][2] = v.z;
	m.data[2][3] = 0;
	m.data[3][0] = 0;
	m.data[3][1] = 0;
	m.data[3][2] = 0;
	m.data[3][3] = 1;
	return m;
}

Matrix4 rotate_matrix_x(float angle)
{
	Matrix4 m;
	m.data[0][0] = 1;
	m.data[1][0] = 0;
	m.data[2][0] = 0;
	m.data[3][0] = 0;
	m.data[0][1] = 0;
	m.data[1][1] = cos(angle);
	m.data[2][1] = -sin(angle);
	m.data[3][1] = 0;
	m.data[0][2] = 0;
	m.data[1][2] = sin(angle);
	m.data[2][2] = cos(angle);
	m.data[3][2] = 0;
	m.data[0][3] = 0;
	m.data[1][3] = 0;
	m.data[2][3] = 0;
	m.data[3][3] = 1;
	return m;
}

Matrix4 rotate_matrix_y(float angle)
{
	Matrix4 m;
	m.data[0][0] = cos(angle);
	m.data[1][0] = 0;
	m.data[2][0] = sin(angle);
	m.data[3][0] = 0;
	m.data[0][1] = 0;
	m.data[1][1] = 1;
	m.data[2][1] = 0;
	m.data[3][1] = 0;
	m.data[0][2] = -sin(angle);
	m.data[1][2] = 0;
	m.data[2][2] = cos(angle);
	m.data[3][2] = 0;
	m.data[0][3] = 0;
	m.data[1][3] = 0;
	m.data[2][3] = 0;
	m.data[3][3] = 1;
	return m;
}

Matrix4 rotate_matrix_z(float angle)
{
	Matrix4 m;
	m.data[0][0] = cos(angle);
	m.data[1][0] = -sin(angle);
	m.data[2][0] = 0;
	m.data[3][0] = 0;
	m.data[0][1] = sin(angle);
	m.data[1][1] = cos(angle);
	m.data[2][1] = 0;
	m.data[3][1] = 0;
	m.data[0][2] = 0;
	m.data[1][2] = 0;
	m.data[2][2] = 1;
	m.data[3][2] = 0;
	m.data[0][3] = 0;
	m.data[1][3] = 0;
	m.data[2][3] = 0;
	m.data[3][3] = 1;
	return m;
}

Matrix4 transpose(Matrix4 m)
{
	Matrix4 r;
	r.data[0][0] = m.data[0][0];
	r.data[0][1] = m.data[1][0];
	r.data[0][2] = m.data[2][0];
	r.data[0][3] = m.data[3][0];
	r.data[1][0] = m.data[0][1];
	r.data[1][1] = m.data[1][1];
	r.data[1][2] = m.data[2][1];
	r.data[1][3] = m.data[3][1];
	r.data[2][0] = m.data[0][2];
	r.data[2][1] = m.data[1][2];
	r.data[2][2] = m.data[2][2];
	r.data[2][3] = m.data[3][2];
	r.data[3][0] = m.data[0][3];
	r.data[3][1] = m.data[1][3];
	r.data[3][2] = m.data[2][3];
	r.data[3][3] = m.data[3][3];
	return r;
}

Matrix4 identity_matrix(void)
{
	Matrix4 m;
	m.data[0][0] = 1;
	m.data[0][1] = 0;
	m.data[0][2] = 0;
	m.data[0][3] = 0;
	m.data[1][0] = 0;
	m.data[1][1] = 1;
	m.data[1][2] = 0;
	m.data[1][3] = 0;
	m.data[2][0] = 0;
	m.data[2][1] = 0;
	m.data[2][2] = 1;
	m.data[2][3] = 0;
	m.data[3][0] = 0;
	m.data[3][1] = 0;
	m.data[3][2] = 0;
	m.data[3][3] = 1;
	return m;
}

Vector4 ldotv(Vector4 v, Matrix4 m)
{
	Vector4 r;
	r.x = m.data[0][0] * v.x + m.data[0][1] * v.y + m.data[0][2] * v.z + m.data[0][3] * v.w;
	r.y = m.data[1][0] * v.x + m.data[1][1] * v.y + m.data[1][2] * v.z + m.data[1][3] * v.w;
	r.z = m.data[2][0] * v.x + m.data[2][1] * v.y + m.data[2][2] * v.z + m.data[2][3] * v.w;
	r.w = m.data[3][0] * v.x + m.data[3][1] * v.y + m.data[3][2] * v.z + m.data[3][3] * v.w;
	return r;
}

Vector4 rdotv(Matrix4 m, Vector4 v)
{
	Vector4 r;
	r.x = m.data[0][0] * v.x + m.data[1][0] * v.y + m.data[2][0] * v.z + m.data[3][0] * v.w;
	r.y = m.data[0][1] * v.x + m.data[1][1] * v.y + m.data[2][1] * v.z + m.data[3][1] * v.w;
	r.z = m.data[0][2] * v.x + m.data[1][2] * v.y + m.data[2][2] * v.z + m.data[3][2] * v.w;
	r.w = m.data[0][3] * v.x + m.data[1][3] * v.y + m.data[2][3] * v.z + m.data[3][3] * v.w;
	return r;
}

Matrix4 dotm(Matrix4 a, Matrix4 b)
{
	Matrix4 r;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			r.data[i][j] = 0;
			for (int k = 0; k < 4; k++)
				r.data[i][j] += a.data[k][j] * b.data[i][k];
		}
	return r;
}

float dotv(Vector3 u, Vector3 v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vector3 mulv(Vector3 a, Vector3 b)
{
	return (Vector3) {
		.x = a.x * b.x,
		.y = a.y * b.y,
		.z = a.z * b.z,
	};
}

Matrix4 lookat_matrix(Vector3 eye, Vector3 center, Vector3 up)
{
	Vector3 forward = combine(center, eye, 1, -1);
	forward = normalize(forward);

	Vector3 right = cross(forward, up);
	right = normalize(right);

	up = cross(right, forward);
	up = normalize(up);

	Matrix4 m;
	m.data[0][0] = right.x;
	m.data[0][1] = up.x;
	m.data[0][2] = -forward.x;
	m.data[0][3] = 0;
	m.data[1][0] = right.y;
	m.data[1][1] = up.y;
	m.data[1][2] = -forward.y;
	m.data[1][3] = 0;
	m.data[2][0] = right.z;
	m.data[2][1] = up.z;
	m.data[2][2] = -forward.z;
	m.data[2][3] = 0;
	m.data[3][0] = -dotv(right, eye);
	m.data[3][1] = -dotv(up, eye);
	m.data[3][2] = dotv(forward, eye);
	m.data[3][3] = 1;
	return m;
}

Matrix4 ortho_matrix(float left, float right, float bottom, float top, float near, float far)
{
	Matrix4 m;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			m.data[i][j] = 0;

	m.data[0][0] = (float)  2 / (right - left);
	m.data[1][1] = (float)  2 / (top - bottom);
	m.data[2][2] = (float) -2 / (far - near);
	m.data[3][3] = 1;

	m.data[3][0] = -(right + left) / (right - left);
	m.data[3][1] = -(top + bottom) / (top - bottom);
	m.data[3][2] = -(far + near) / (far - near);

	return m;
}

Matrix4 perspective_matrix(float fov, float aspect, float near, float far)
{
	Matrix4 m;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			m.data[i][j] = 0;

	float tan_half_fov = tanf(fov / 2.0);

	m.data[0][0] = (float) 1 / (aspect * tan_half_fov);
	m.data[1][1] = (float) 1 / tan_half_fov;
	m.data[2][2] = -(far + near) / (far - near);
	m.data[2][3] = -1;
	m.data[3][2] = -(2 * far * near) / (far - near);

	return m;
}

void print_matrix(Matrix4 m)
{
	printf(
		"| %f %f %f %f |\n"
		"| %f %f %f %f |\n"
		"| %f %f %f %f |\n"
		"| %f %f %f %f |\n\n",
		m.data[0][0],
		m.data[0][1],
		m.data[0][2],
		m.data[0][3],
		m.data[1][0],
		m.data[1][1],
		m.data[1][2],
		m.data[1][3],
		m.data[2][0],
		m.data[2][1],
		m.data[2][2],
		m.data[2][3],
		m.data[3][0],
		m.data[3][1],
		m.data[3][2],
		m.data[3][3]);
}

void print_matrix3(Matrix3 m)
{
	printf(
		"| %f %f %f |\n"
		"| %f %f %f |\n"
		"| %f %f %f |\n\n",
		m.data[0][0],
		m.data[0][1],
		m.data[0][2],
		m.data[1][0],
		m.data[1][1],
		m.data[1][2],
		m.data[2][0],
		m.data[2][1],
		m.data[2][2]);
}

static bool gluInvertMatrix(const float m[16], float invOut[16])
{
    float inv[16], det;
    int i;

    inv[0] =
		m[5]  * m[10] * m[15] - 
		m[5]  * m[11] * m[14] - 
		m[9]  * m[6]  * m[15] + 
		m[9]  * m[7]  * m[14] +
		m[13] * m[6]  * m[11] - 
		m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}

bool invert(Matrix4 a, Matrix4 *b)
{
	return gluInvertMatrix((float*) &a, (float*) b);
}
