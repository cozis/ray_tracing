#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "mesh.h"

#define TINYOBJ_LOADER_C_IMPLEMENTATION
#include "tinyobj_loader_c.h"

void append_vertex(VertexArray *array, Vertex v)
{
	if (array->size == array->capacity) {
		if (array->capacity == 0) {
			array->data = malloc(8 * sizeof(Vertex));
			array->capacity = 8;
		} else {
			array->data = realloc(array->data, 2 * array->capacity * sizeof(Vertex));
			array->capacity *= 2;
		}
		if (!array->data) {
			printf("OUT OF MEMORY\n");
			abort();
		}
	}
	array->data[array->size++] = v;
}

static Vector3 get_sphere_point(float angle_x, float angle_y, float radius)
{
    Vector3 p;
    p.x = radius * sin(angle_y) * cos(2 * angle_x);
    p.y = radius * cos(angle_y);
    p.z = radius * sin(angle_y) * sin(2 * angle_x);
    return p;
}

typedef struct {
    Vertex a;
    Vertex b;
    Vertex c;
} Triangle;

static void calculate_and_set_normals(Triangle *T)
{
    #define X1 (T->b.x - T->a.x)
    #define Y1 (T->b.y - T->a.y)
    #define Z1 (T->b.z - T->a.z)

    #define X2 (T->c.x - T->a.x)
    #define Y2 (T->c.y - T->a.y)
    #define Z2 (T->c.z - T->a.z)

    /*
     * xnormal = y1*z2 - z1*y2
     * ynormal = z1*x2 - x1*z2
     * znormal = x1*y2 - y1*x2
    */

    Vector3 n;
    n.x = Y1 * Z2 - Z1 * Y2;
    n.y = Z1 * X2 - X1 * Z2;
    n.z = X1 * Y2 - Y1 * X2;

    T->a.nx = n.x;
    T->a.ny = n.y;
    T->a.nz = n.z;

    T->b.nx = n.x;
    T->b.ny = n.y;
    T->b.nz = n.z;

    T->c.nx = n.x;
    T->c.ny = n.y;
    T->c.nz = n.z;

    #undef X1
    #undef Y1
    #undef Z1
    #undef X2
    #undef Y2
    #undef Z2
}

static Triangle make_triangle(Vector3 a, Vector3 b, Vector3 c)
{
    Triangle T;
    T.a = (Vertex) {a.x, a.y, a.z};
    T.b = (Vertex) {b.x, b.y, b.z};
    T.c = (Vertex) {c.x, c.y, c.z};
    calculate_and_set_normals(&T);
    return T;
}

VertexArray make_sphere_mesh_2(float radius, int num_segms, bool fake_normals)
{
    VertexArray vertices = {0, 0, 0};

    int x_num_segms = num_segms;
    int y_num_segms = num_segms;
    for (int i = 0; i < x_num_segms; i++)
        for (int j = 0; j < y_num_segms; j++) {

            int g = j;
            if (g == y_num_segms-1)
                g = 0;

            Vector3 p1 = get_sphere_point((i + 0) * 3.14 / x_num_segms, (g + 0) * 3.14 / y_num_segms, radius);
            Vector3 p2 = get_sphere_point((i + 1) * 3.14 / x_num_segms, (g + 0) * 3.14 / y_num_segms, radius);
            Vector3 p3 = get_sphere_point((i + 1) * 3.14 / x_num_segms, (g + 1) * 3.14 / y_num_segms, radius);
            Vector3 p4 = get_sphere_point((i + 0) * 3.14 / x_num_segms, (g + 1) * 3.14 / y_num_segms, radius);

            Triangle t1 = make_triangle(p1, p2, p3);

            if (fake_normals) {
                t1.a.nx = p1.x;
                t1.a.ny = p1.y;
                t1.a.nz = p1.z;

                t1.b.nx = p2.x;
                t1.b.ny = p2.y;
                t1.b.nz = p2.z;

                t1.c.nx = p3.x;
                t1.c.ny = p3.y;
                t1.c.nz = p3.z;
            }

			append_vertex(&vertices, t1.a);
			append_vertex(&vertices, t1.b);
			append_vertex(&vertices, t1.c);

            Triangle t2 = make_triangle(p4, p1, p3);

            if (fake_normals) {
                t2.a.nx = p4.x;
                t2.a.ny = p4.y;
                t2.a.nz = p4.z;

                t2.b.nx = p1.x;
                t2.b.ny = p1.y;
                t2.b.nz = p1.z;

                t2.c.nx = p3.x;
                t2.c.ny = p3.y;
                t2.c.nz = p3.z;
            }

			append_vertex(&vertices, t2.a);
			append_vertex(&vertices, t2.b);
			append_vertex(&vertices, t2.c);
        }
    return vertices;
}

VertexArray make_sphere_mesh(float radius)
{
	return make_sphere_mesh_2(radius, 32, true);
}

VertexArray make_cube_mesh(void)
{
    float vertices[] = {

        1.0, 1.0, 0.0,   0.0, 0.0, -1.0,
        1.0, 0.0, 0.0,   0.0, 0.0, -1.0,
        0.0, 0.0, 0.0,   0.0, 0.0, -1.0,
        0.0, 1.0, 0.0,   0.0, 0.0, -1.0,
        1.0, 1.0, 0.0,   0.0, 0.0, -1.0,
        0.0, 0.0, 0.0,   0.0, 0.0, -1.0,

        0.0, 0.0, 1.0,   0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,   0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,   0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,   0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,   0.0, 0.0, 1.0,
        0.0, 1.0, 1.0,   0.0, 0.0, 1.0,

        0.0, 0.0, 0.0,   0.0, -1.0, 0.0,
        1.0, 0.0, 0.0,   0.0, -1.0, 0.0,
        1.0, 0.0, 1.0,   0.0, -1.0, 0.0,
        0.0, 0.0, 0.0,   0.0, -1.0, 0.0,
        1.0, 0.0, 1.0,   0.0, -1.0, 0.0,
        0.0, 0.0, 1.0,   0.0, -1.0, 0.0,

        1.0, 1.0, 1.0,   0.0, 1.0, 0.0,
        1.0, 1.0, 0.0,   0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,   0.0, 1.0, 0.0,
        0.0, 1.0, 1.0,   0.0, 1.0, 0.0,
        1.0, 1.0, 1.0,   0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,   0.0, 1.0, 0.0,

        0.0, 1.0, 1.0,  -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,  -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,  -1.0, 0.0, 0.0,

        0.0, 0.0, 1.0,  -1.0, 0.0, 0.0,
        0.0, 1.0, 1.0,  -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,  -1.0, 0.0, 0.0,

        1.0, 0.0, 0.0,  1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,  1.0, 0.0, 0.0,
        1.0, 1.0, 1.0,  1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,  1.0, 0.0, 0.0,
        1.0, 1.0, 1.0,  1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,  1.0, 0.0, 0.0,
    };

   VertexArray result = {0, 0, 0}; 
   for (int i = 0; i < (int) (sizeof(vertices)/sizeof(vertices[0])); i += 6) {
        Vertex v;
        v.x = vertices[i + 0];
        v.y = vertices[i + 1];
        v.z = vertices[i + 2];
        v.nx = vertices[i + 3];
        v.ny = vertices[i + 4];
        v.nz = vertices[i + 5];
        v.tx = 0;
        v.ty = 0;
		append_vertex(&result, v);
   }

   return result;
}

static void
get_file_data_callback(void *context, const char *filename, 
	const int is_mtl, const char *obj_filename, char **data, size_t *len)
{
	(void) context;

	if (filename == NULL) {
		fprintf(stderr, "null filename\n");
		(*data) = NULL;
		(*len) = 0;
		return;
	}

	size_t data_len = 0;

	printf("Reading file '%s'\n", filename);
	*data = load_file(filename, &data_len);
	(*len) = data_len;

	char **free_me = context;
	if (*free_me == NULL) *free_me = *data;
}

bool load_mesh_from_file(const char *file, VertexArray *result)
{

	tinyobj_attrib_t attrib;

	tinyobj_shape_t* shapes = NULL;
	size_t num_shapes;

	tinyobj_material_t* materials = NULL;
	size_t num_materials;

	char *free_me = NULL;

	unsigned int flags = TINYOBJ_FLAG_TRIANGULATE;
	int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &materials, &num_materials, file, get_file_data_callback, &free_me, flags);
	if (ret != TINYOBJ_SUCCESS) {
		if (free_me) free(free_me);
		printf("Failed loading '%s'\n", file);
		return false;
	}

	*result = (VertexArray) {0, 0, 0};

	size_t face_offset = 0;
	for (int i = 0; i < (int) attrib.num_face_num_verts; i++) {

		assert(attrib.face_num_verts[i] % 3 == 0); /* assume all triangle faces. */
		for (size_t f = 0; f < (size_t)attrib.face_num_verts[i] / 3; f++) {

			tinyobj_vertex_index_t idx0 = attrib.faces[face_offset + f * 3 + 0];
			tinyobj_vertex_index_t idx1 = attrib.faces[face_offset + f * 3 + 1];
			tinyobj_vertex_index_t idx2 = attrib.faces[face_offset + f * 3 + 2];

			Vertex v0;
			Vertex v1;
			Vertex v2;

			/*
			 * Positions
			 */

			v0.x = attrib.vertices[idx0.v_idx * 3 + 0];
			v0.y = attrib.vertices[idx0.v_idx * 3 + 1];
			v0.z = attrib.vertices[idx0.v_idx * 3 + 2];

			v1.x = attrib.vertices[idx1.v_idx * 3 + 0];
			v1.y = attrib.vertices[idx1.v_idx * 3 + 1];
			v1.z = attrib.vertices[idx1.v_idx * 3 + 2];

			v2.x = attrib.vertices[idx2.v_idx * 3 + 0];
			v2.y = attrib.vertices[idx2.v_idx * 3 + 1];
			v2.z = attrib.vertices[idx2.v_idx * 3 + 2];

			/*
			 * Normals
			 */

			v0.nx = attrib.normals[idx0.vn_idx * 3 + 0];
			v0.ny = attrib.normals[idx0.vn_idx * 3 + 1];
			v0.nz = attrib.normals[idx0.vn_idx * 3 + 2];

			v1.nx = attrib.normals[idx1.vn_idx * 3 + 0];
			v1.ny = attrib.normals[idx1.vn_idx * 3 + 1];
			v1.nz = attrib.normals[idx1.vn_idx * 3 + 2];

			v2.nx = attrib.normals[idx2.vn_idx * 3 + 0];
			v2.ny = attrib.normals[idx2.vn_idx * 3 + 1];
			v2.nz = attrib.normals[idx2.vn_idx * 3 + 2];

			/*
			 * Texture coordinates
			 */

			v0.tx = attrib.normals[idx0.vt_idx * 2 + 0];
			v0.ty = attrib.normals[idx0.vt_idx * 2 + 1];

			v1.tx = attrib.normals[idx1.vt_idx * 2 + 0];
			v1.ty = attrib.normals[idx1.vt_idx * 2 + 1];

			v2.tx = attrib.normals[idx2.vt_idx * 2 + 0];
			v2.ty = attrib.normals[idx2.vt_idx * 2 + 1];

			append_vertex(result, v0);
			append_vertex(result, v1);
			append_vertex(result, v2);
		}

		face_offset += (size_t)attrib.face_num_verts[i];
	}

	if (free_me)
		free(free_me);

	tinyobj_attrib_free(&attrib);
	tinyobj_shapes_free(shapes, num_shapes);
	tinyobj_materials_free(materials, num_materials);
	return true;
}