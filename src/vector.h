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
#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

#include <stdbool.h>

typedef struct {
	float x;
	float y;
	float z;
} Vector3;

typedef struct {
	float x;
	float y;
	float z;
	float w;
} Vector4;

typedef struct {
	float data[4][4];
} Matrix4;

typedef struct {
	float data[3][3];
} Matrix3;

typedef struct {
	Vector3 origin;
	Vector3 direction;
} Ray;

typedef struct {
	Vector3 center;
	float   radius;
} Sphere;

float deg2rad(float deg);

float maxf(float x, float y);
float minf(float x, float y);
float absf(float x);
float clamp(float x, float min, float max);
Vector3 maxv(Vector3 a, Vector3 b);
bool isnanv(Vector3 v);
bool iszerof(float f);
bool iszerov(Vector3 v);
float avgv(Vector3 v);

Vector3 vec_from_scalar(float s);

Matrix4 translate_matrix(Vector3 v, float f);
Matrix4 identity_matrix(void);
Matrix4 scale_matrix(Vector3 v);
Matrix4 rotate_matrix_x(float angle);
Matrix4 rotate_matrix_y(float angle);
Matrix4 rotate_matrix_z(float angle);
Matrix4 lookat_matrix(Vector3 pos, Vector3 front, Vector3 up);
Matrix4 ortho_matrix(float left, float right, float bottom, float top, float near, float far);
Matrix4 perspective_matrix(float fov, float aspect, float near, float far);

void print_matrix(Matrix4 m);

float   norm2_of(Vector3 v);
float   norm_of(Vector3 v);
Vector3 normalize(Vector3 v);
Vector3 scalev(Vector3 v, float f);
Vector3 combine(Vector3 u, Vector3 v, float a, float b);
Vector3 combine4(Vector3 u, Vector3 v, Vector3 g, Vector3 t, float a, float b, float c, float d);
Vector3 cross(Vector3 u, Vector3 v);

Vector3 mulv(Vector3 a, Vector3 b);

Vector4 ldotv(Vector4 v, Matrix4 m);
Vector4 rdotv(Matrix4 m, Vector4 v);
float   dotv(Vector3 u, Vector3 v);
Matrix4 dotm(Matrix4 a, Matrix4 b);
Matrix4 transpose(Matrix4 m);
bool invert(Matrix4 a, Matrix4 *inv);

Vector3 random_vector(void);
Vector3 random_direction(void);
Vector3 reflect(Vector3 dir, Vector3 normal);

#ifndef M_PI
#define M_PI 3.1415926538
#endif

#endif