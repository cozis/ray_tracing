#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "vector.h"

#define EPSILON 0.00001

bool isnanv(Vector3 v)
{
	return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

float deg2rad(float deg)
{
	return 3.14159265358979323846 * deg / 180;
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

/*
void lina_transpose(float *A, float *B, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    if(m == 1 || n == 1) {
        // For a matrix with height or width of 1
        // row-major and column-major order coincide,
        // so the stransposition doesn't change the
        // the memory representation. A simple copy
        // does the job.

            if(A != B) // Does the copy or the branch cost more?
                memcpy(B, A, sizeof(A[0]) * m * n);
    
    } else if(m == n) {

        // Iterate over the upper triangular portion of
        // the matrix and switch each element with the
        // corresponding one in the lower triangular portion.
        // NOTE: We're assuming A,B might be the same matrix.
        //       If A,B are the same matrix, then the diagonal
        //       is copied onto itself. By removing the +1 in
        //       the inner loop, the copying of the diagonal
        //       is avoided.

        for(int i = 0; i < n; i += 1)
            for(int j = 0; j < i+1; j += 1) {
                float temp = A[i*n + j];
                B[i*n + j] = A[j*n + i];
                B[j*n + i] = temp;
            }

    } else {
        // Not only the matrix needs to be transposed
        // assuming the destination matrix is the same
        // as the source matrix, but the memory representation
        // of the matrix needs to switch from row-major
        // to col-major, so it's not as simple as switching
        // value's positions.
        // This algorithm starts from the A[0][1] value and
        // moves it where it needs to go, then gets the value
        // that was at that position and puts that in it's
        // new position. This process is iterated until the
        // starting point A[0][1] is overwritten with the
        // new value. In this process the first and last
        // value of the matrix never move.

        B[0] = A[0];
        B[m*n - 1] = A[m*n - 1];

        float item = A[1];
        int    next = m;

        while(next != 1) {
            float temp = A[next];
            B[next] = item;
            item = temp;
            next = (next % n) * m + (next / n);
        }

        B[1] = item;
    }
}

int lina_decompLUP(float *A, float *L, 
                   float *U, int   *P, 
                   int n)
{
    assert(n > 0);
    assert(A != L && A != U && L != U);

    for (int i = 0; i < n; i++)
        P[i] = i;

    int swaps = 0;
    for (int i = 0; i < n; i++) {

        int v = P[i];
        float max_v = A[v * n + i];
        int    max_i = i;
        
        for (int j = i+1; j < n; j++) {
            int u = P[j];
            float abs = fabs(A[u * n + j]);
            if (abs > max_v) {
                max_v = abs;
                max_i = j;
            }
        }

        if (max_i != i) {

            // Swap rows
            int temp = P[i];
            P[i] = P[max_i];
            P[max_i] = temp;

            swaps++;
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U[i * n + j] = A[P[i] * n + j];

    memset(L, 0, sizeof(float) * n * n);
    for (int i = 0; i < n; i++)
        L[i * n + i] = 1;

    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            float u = U[i * n + i];
            L[j * n + i] = U[j * n + i] / u;
            for (int k = 0; k < n; k++)
                U[j * n + k] -= L[j * n + i] * U[i * n + k];
        }

    return swaps;
}

// Function: lina_det
//
//   Calculates the determinant of the n by n matrix A
//   and returns it throught the output parameter [det].
//
//   If not enough memory is available, false is returned,
//   else true is returned.
//
// Notes:
//   - The output parameter [det] is optional. (you can
//     ignore the result by passing NULL).
//
bool lina_det(float *A, int n, float *det)
{
    // Allocate the space for the L,U matrices.
    // I can't think of a version of this algorithm
    // where a temporary buffer isn't necessary.
    float *T = (float*) malloc(sizeof(float) * n * n * 2 + sizeof(int) * n);
    if (T == NULL)
        return false;

    // Do the decomposition
    float *L = T;
    float *U = L + (n * n);
    int    *P = (int*) (U + (n * n));
    
    int swaps = lina_decompLUP(A, L, U, P, n);
    if (swaps < 0) {
        free(T);
        return false;
    }

    // Knowing that
    //
    //   A = LU
    //
    // then
    //
    //   det(A) = det(LU) = det(L)det(U)
    //
    // Since L and U are triangular, their 
    // determinant is the product of their 
    // diagonals, so the product of the 
    // determinants is the product of both 
    // the diagonals.

    float prod = 1;
    for (int i = 0; i < n; i++) {
        float l = L[i * n + i];
        float u = U[i * n + i];
        prod *= l * u;
    }

    if (swaps & 1)
        prod = -prod;

    if (det)
        *det = prod;

    free(T);
    return true;
}

// Create the n-1 by n-1 matrix D obtained by
// removing the [del_col] column and [del_row]
// frow the n by n matrix M.
static void 
copyMatrixWithoutRowAndCol(float *M, float *D, int n, 
                           int del_col, int del_row)
{
    // Copy the upper-left portion of matrix M
    // that comes before the deleted column and
    // row.
    for (int i = 0; i < del_row; i++)
        for (int j = 0; j < del_col; j++)
            D[i * (n-1) + j] = M[i * n + j];

    // Copy the lower left portion that comes
    // after both the deleted column and row.
    for (int i = del_row+1; i < n; i++)
        for (int j = del_col+1; j < n; j++)
            D[(i-1) * (n-1) + (j-1)] = M[i * n + j];

    // Copy the bottom portion that comes after
    // the deleted row but before the deleted column.
    for (int i = del_row+1; i < n; i++)
        for (int j = 0; j < del_col; j++)
            D[(i-1) * (n-1) + j] = M[i * n + j];

    // Copy the right portion that comes after
    // the deleted column but before the deleted row.
    for (int i = 0; i < del_row; i++)
        for (int j = del_col+1; j < n; j++)
            D[i * (n-1) + (j-1)] = M[i * n + j];
}

bool lina_inverse(Matrix4 M, Matrix4 *D)
{
    float det;
    if (!lina_det((float*) &M, 4, &det))
        return false;

	printf("det=%f\n", det);

    if (det == 0)
        return false; // The matrix can't be inverted

	Matrix3 T;
    Matrix4 M_t = transpose(M);

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            
            copyMatrixWithoutRowAndCol((float*) &M_t, (float*)&T, 4, j, i);

            float det2;
            if (!lina_det((float*) &T, 3, &det2))
                return false;

			printf("-----------------\n");

			print_matrix3(T);
			printf("det2=%f\n", det2);

            // If the determinant of M isn't zero,
            // neither is this!
            assert(det2 != 0);

            bool i_is_odd = i & 1;
            bool j_is_odd = j & 1;
            int sign = (i_is_odd == j_is_odd) ? 1 : -1;

            D->data[i][j] = sign * det2 / det;
        }
    return true;
}

Matrix4 invert(Matrix4 m)
{
	Matrix4 r;
	if (!lina_inverse(m, &r)) {
		printf("Couldn't invert!\n");
		abort();
	}
	return r;
}
*/

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
