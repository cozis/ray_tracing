#include "vector.h"

#define MAX_OBJECTS 1024

typedef struct {
	Vector3 albedo;
	float   roughness;
	float   reflectance;
	float   metallic;
	float   emission_power;
	Vector3 emission_color;
} Material;

typedef struct {
	Vector3 origin;
	Vector3 size;
} Cube;

typedef enum {
	OBJECT_CUBE,
	OBJECT_SPHERE,
} ObjectType;

typedef struct {
	ObjectType type;
	union {
		Sphere sphere;
		Cube cube;
	};
	Material material;
} Object;

typedef struct {
	Object objects[MAX_OBJECTS];
	int num_objects;
} Scene;

typedef struct {
	float   distance;
	Vector3 point;
	Vector3 normal;
	int     object;
} HitInfo;

Vector3 origin_of(Object o);
HitInfo trace_ray(Ray ray, Scene *scene);
bool    parse_scene_file(char *file, Scene *scene);