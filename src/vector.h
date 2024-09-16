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
Vector3 scale(Vector3 v, float f);
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

#endif