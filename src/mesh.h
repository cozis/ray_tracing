#include "vector.h"

typedef struct {
    float x, y, z;
    float nx, ny, nz;
    float tx, ty;
} Vertex;

typedef struct {
	Vertex *data;
	int size;
	int capacity;
} VertexArray;

void append_vertex(VertexArray *array, Vertex v);

VertexArray make_sphere_mesh(float radius);
VertexArray make_sphere_mesh_2(float radius, int num_segms, bool fake_normals);
VertexArray make_cube_mesh(void);

bool load_mesh_from_file(const char *file, VertexArray *result);