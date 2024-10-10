#include <math.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "utils.h"
#include "scene.h"

Vector3 origin_of(Object o)
{
	if (o.type == OBJECT_SPHERE)
		return o.sphere.center;
	return combine(o.cube.origin, o.cube.size, 1, 0.5);
}

static bool intersect_cube(Ray r, Cube c, float *tnear, float *tfar, Vector3 *normal)
{
	float txmin, txmax;
	float tymin, tymax;
	float tzmin, tzmax;

	float tn;
	float tf;

	Vector3 a = c.origin;
	Vector3 b = combine(c.origin, c.size, 1, 1);

	int hit_axis = 0; // 0=x, 1=y, 2=z

	if (r.direction.x >= 0) {
		txmin = (a.x - r.origin.x) / r.direction.x;
		txmax = (b.x - r.origin.x) / r.direction.x;
	} else {
		txmax = (a.x - r.origin.x) / r.direction.x;
		txmin = (b.x - r.origin.x) / r.direction.x;
	}

	if (r.direction.y >= 0) {
		tymin = (a.y - r.origin.y) / r.direction.y;
		tymax = (b.y - r.origin.y) / r.direction.y;
	} else {
		tymax = (a.y - r.origin.y) / r.direction.y;
		tymin = (b.y - r.origin.y) / r.direction.y;
	}

	if (txmin > tymax || tymin > txmax)
		return false;

	if (tymin > txmin) { txmin = tymin; hit_axis = 1; }
	if (tymax < txmax) txmax = tymax;

	if (r.direction.z >= 0) {
		tzmin = (a.z - r.origin.z) / r.direction.z;
		tzmax = (b.z - r.origin.z) / r.direction.z;
	} else {
		tzmax = (a.z - r.origin.z) / r.direction.z;
		tzmin = (b.z - r.origin.z) / r.direction.z;
	}

	if (txmin > tzmax || tzmin > txmax)
		return false;
	
	if (tzmin > txmin) { txmin = tzmin; hit_axis = 2; };
	if (tzmax < txmax) txmax = tzmax;

	if (tnear) *tnear = txmin;
	if (tfar)  *tfar  = txmax;
	if (normal) {
		switch (hit_axis) {
			case 0: *normal = r.direction.x > 0 ? (Vector3) {-1, 0, 0} : (Vector3) {1, 0, 0}; break;
			case 1: *normal = r.direction.y > 0 ? (Vector3) {0, -1, 0} : (Vector3) {0, 1, 0}; break;
			case 2: *normal = r.direction.z > 0 ? (Vector3) {0, 0, -1} : (Vector3) {0, 0, 1}; break;
		}
	}
	return true;
}

static bool intersect_sphere(Ray r, Sphere s, float *t)
{
	/*
	 * Any point of the ray can be written as
	 *
	 *     P(t) = O + t * D
	 * 
	 * with O origin and D direction.
	 * 
	 * All points P=(x,y,z) of a sphere can be described as
	 * those (and only those) that satisfy the equation
	 *
	 *     x^2 + y^2 + z^2 = R^2
	 *     P^2 - R^2 = 0
	 * 
	 * with R radius of the sphere. The sphere here is centered
	 * at the origin.
	 * 
	 * Intersection points of the ray with the sphere must satisfy
	 * both:
	 * 
	 *     P(t) = O + t * D
	 *     P^2 - R^2 = 0
	 * 
	 *     => (O + tD)^2 - R^2 = 0
	 *     => t^2 * D^2 + t * 2OD + O^2 - R^2 = 0
	 * 
	 * we can use the quadratic formula here, and more specifically
	 * the discriminant to check if solutions exist and how many
	 */
	Vector3 oc = combine(s.center, r.origin, 1, -1);
	float a = dotv(r.direction, r.direction);
	float b = -2 * dotv(oc, r.direction);
	float c = dotv(oc, oc) - s.radius * s.radius;

	float discr = b*b - 4*a*c;

	if (discr > 0) {
		float s0 = (- b + sqrt(discr)) / (2 * a);
		float s1 = (- b - sqrt(discr)) / (2 * a);
		if (s0 > s1) {
			float tmp = s0;
			s0 = s1;
			s1 = tmp;
		}
		if (s0 < 0) {
			s0 = s1;
			if (s0 < 0) return false;
		}
		if (t) *t = s0;
		return true;
	}

	// Zero solutions
	return false;
}

static bool intersect_object(Ray r, Object o, float *t, Vector3 *normal)
{
	switch (o.type) {

		case OBJECT_CUBE:
		return intersect_cube(r, o.cube, t, NULL, normal);

		case OBJECT_SPHERE:
		if (intersect_sphere(r, o.sphere, t)) {
			if (normal) {
				Vector3 hit_point = combine(r.origin, r.direction, 1, *t);
				*normal = normalize(combine(hit_point, o.sphere.center, 1, -1));
			}
			return true;
		}
		return false;
	}
	return false;
}

HitInfo trace_ray(Ray ray, Scene *scene)
{
	ray.direction = normalize(ray.direction);

	float   nearest_t = FLT_MAX;
	int     nearest_object = -1;
	Vector3 nearest_normal;
	for (int i = 0; i < scene->num_objects; i++) {
		float t;
		Vector3 n;
		if (!intersect_object(ray, scene->objects[i], &t, &n))
			continue;
		if (t >= 0 && t < nearest_t) {
			nearest_t = t;
			nearest_object = i;
			nearest_normal = n;
		}
	}

	if (nearest_object == -1) {
		HitInfo result;
		result.distance = -1;
		result.normal   = (Vector3) {0, 0, 0};
		result.point    = (Vector3) {0, 0, 0};
		result.object   = -1;
		return result;
	} else {
		HitInfo result;
		result.distance = nearest_t;
		result.normal   = nearest_normal;
		result.point    = combine(ray.origin, ray.direction, 1, nearest_t);
		result.object   = nearest_object;
		return result;
	}
}


typedef enum {
	PROP_ALBEDO,
	PROP_ROUGHNESS,
	PROP_REFLECTANCE,
	PROP_METALLIC,
	PROP_EMISSION_POWER,
	PROP_EMISSION_COLOR,
	PROP_RADIUS,
	PROP_CENTER,
	PROP_ORIGIN,
	PROP_SIZE,
} Property;

static bool parse_scene_string(char *src, size_t len, Scene *scene)
{
	scene->num_objects = 0;

	int line = 1;
	size_t i = 0;
	for (;;) {

		while (i < len && is_space(src[i])) {
			if (src[i] == '\n') line++;
			i++;
		}

		if (i == len)
			break;
		
		Object object;

		if (5 < len - i
			&& src[i+0] == 's'
			&& src[i+1] == 'p'
			&& src[i+2] == 'h'
			&& src[i+3] == 'e'
			&& src[i+4] == 'r'
			&& src[i+5] == 'e') {
			object.type = OBJECT_SPHERE;
			object.sphere.center = (Vector3) {0, 0, 0};
			object.sphere.radius = 1;
			object.material.albedo = (Vector3) {0.44, 0.68, 0.84};
			object.material.roughness = 0;
			object.material.reflectance = 0.2;
			object.material.metallic = 0;
			object.material.emission_power = 0;
			object.material.emission_color = (Vector3) {1, 1, 1};
			i += 6;
		} else if (3 < len - i
			&& src[i+0] == 'c'
			&& src[i+1] == 'u'
			&& src[i+2] == 'b'
			&& src[i+3] == 'e') {
			object.type = OBJECT_CUBE;
			object.cube.origin = (Vector3) {0, 0, 0};
			object.cube.size = (Vector3) {1, 1, 1};
			object.material.albedo = (Vector3) {0.44, 0.68, 0.84};
			object.material.roughness = 0;
			object.material.reflectance = 0.2;
			object.material.metallic = 0;
			object.material.emission_power = 0;
			object.material.emission_color = (Vector3) {1, 1, 1};
			i += 4;
		} else {
			fprintf(stderr, "Error: Invalid character (line %d)\n", line);
			return false;
		}

		for (;;) {
			// Skip spaces before property
			while (i < len && is_space(src[i])) {
				if (src[i] == '\n') line++;
				i++;
			}
			
			int valuetype; // 0 for float, 1 for color (3 floats)

			Property prop;
			if (6 < len - i
				&& src[i+0] == 'a'
				&& src[i+1] == 'l'
				&& src[i+2] == 'b'
				&& src[i+3] == 'e'
				&& src[i+4] == 'd'
				&& src[i+5] == 'o') {
				valuetype = 1;
				prop = PROP_ALBEDO;
				i += 9;
			} else if (8 < len - i
				&& src[i+0] == 'r'
				&& src[i+1] == 'o'
				&& src[i+2] == 'u'
				&& src[i+3] == 'g'
				&& src[i+4] == 'h'
				&& src[i+5] == 'n'
				&& src[i+6] == 'e'
				&& src[i+7] == 's'
				&& src[i+8] == 's') {
				valuetype = 0;
				prop = PROP_ROUGHNESS;
				i += 9;
			} else if (10 < len - i
				&& src[i+0] == 'r'
				&& src[i+1] == 'e'
				&& src[i+2] == 'f'
				&& src[i+3] == 'l'
				&& src[i+4] == 'e'
				&& src[i+5] == 'c'
				&& src[i+6] == 't'
				&& src[i+7] == 'a'
				&& src[i+8] == 'n'
				&& src[i+9] == 'c'
				&& src[i+10] == 'e') {
				valuetype = 0;
				prop = PROP_REFLECTANCE;
				i += 11;
			} else if (7 < len - i
				&& src[i+0] == 'm'
				&& src[i+1] == 'e'
				&& src[i+2] == 't'
				&& src[i+3] == 'a'
				&& src[i+4] == 'l'
				&& src[i+5] == 'l'
				&& src[i+6] == 'i'
				&& src[i+7] == 'c') {
				valuetype = 0;
				prop = PROP_METALLIC;
				i += 11;
			} else if (13 < len - i
				&& src[i+0] == 'e'
				&& src[i+1] == 'm'
				&& src[i+2] == 'i'
				&& src[i+3] == 's'
				&& src[i+4] == 's'
				&& src[i+5] == 'i'
				&& src[i+6] == 'o'
				&& src[i+7] == 'n'
				&& src[i+8] == '_'
				&& src[i+9] == 'p'
				&& src[i+10] == 'o'
				&& src[i+11] == 'w'
				&& src[i+12] == 'e'
				&& src[i+13] == 'r') {
				valuetype = 0;
				prop = PROP_EMISSION_POWER;
				i += 14;
			} else if (13 < len - i
				&& src[i+0] == 'e'
				&& src[i+1] == 'm'
				&& src[i+2] == 'i'
				&& src[i+3] == 's'
				&& src[i+4] == 's'
				&& src[i+5] == 'i'
				&& src[i+6] == 'o'
				&& src[i+7] == 'n'
				&& src[i+8] == '_'
				&& src[i+9] == 'c'
				&& src[i+10] == 'o'
				&& src[i+11] == 'l'
				&& src[i+12] == 'o'
				&& src[i+13] == 'r') {
				valuetype = 1;
				prop = PROP_EMISSION_COLOR;
				i += 14;
			} else if (5 < len - i
				&& src[i+0] == 'r'
				&& src[i+1] == 'a'
				&& src[i+2] == 'd'
				&& src[i+3] == 'i'
				&& src[i+4] == 'u'
				&& src[i+5] == 's') {
				if (object.type != OBJECT_SPHERE) {
					fprintf(stderr, "Poperty 'radius' only allowed on spheres (line %d)\n", line);
					return false;
				}
				valuetype = 0;
				prop = PROP_RADIUS;
				i += 6;
			} else if (5 < len - i
				&& src[i+0] == 'c'
				&& src[i+1] == 'e'
				&& src[i+2] == 'n'
				&& src[i+3] == 't'
				&& src[i+4] == 'e'
				&& src[i+5] == 'r') {
				if (object.type != OBJECT_SPHERE) {
					fprintf(stderr, "Poperty 'center' only allowed on spheres (line %d)\n", line);
					return false;
				}
				valuetype = 1;
				prop = PROP_CENTER;
				i += 6;
			} else if (5 < len - i
				&& src[i+0] == 'o'
				&& src[i+1] == 'r'
				&& src[i+2] == 'i'
				&& src[i+3] == 'g'
				&& src[i+4] == 'i'
				&& src[i+5] == 'n') {
				if (object.type != OBJECT_CUBE) {
					fprintf(stderr, "Poperty 'origin' only allowed on cubes (line %d)\n", line);
					return false;
				}
				valuetype = 1;
				prop = PROP_ORIGIN;
				i += 6;
			} else if (3 < len - i
				&& src[i+0] == 's'
				&& src[i+1] == 'i'
				&& src[i+2] == 'z'
				&& src[i+3] == 'e') {
				if (object.type != OBJECT_CUBE) {
					fprintf(stderr, "Poperty 'size' only allowed on cubes (line %d)\n", line);
					return false;
				}
				valuetype = 1;
				prop = PROP_SIZE;
				i += 4;
			} else
				// Not a valid property name
				break;

			// Consume spaces before the value
			while (i < len && is_space(src[i])) {
				if (src[i] == '\n') line++;
				i++;
			}
			if (i == len) {
				fprintf(stderr, "Error: Property value is missing (line %d)\n", line);
				return false;
			}

			float   value0;
			Vector3 value1;
			if (valuetype == 0) {
				// Parse a single float
				int sign = 1;
				if (src[i] == '-') {
					sign = -1;
					i++;
					if (i == len || !is_digit(src[i])) {
						fprintf(stderr, "Error: Missing number after minus sign (line %d)\n", line);
						return false;
					}
				} else if (!is_digit(src[i])) {
					fprintf(stderr, "Error: Missing number after property name (line %d)\n", line);
					return false;
				}
				value0 = 0;
				do {
					int d = src[i] - '0';
					value0 = value0 * 10 + d;
					i++;
				} while (i < len && is_digit(src[i]));
				if (i < len && src[i] == '.') {
					i++; // Skip the dot
					if (i == len || !is_digit(src[i])) {
						fprintf(stderr, "Error: Missing decimal part after dot (line %d)\n", line);
						return false;
					}
					float q = 1.0f / 10;
					do {
						int d = src[i] - '0';
						value0 += q * d;
						q /= 10;
						i++;
					} while (i < len && is_digit(src[i]));
				}
				value0 *= sign;
			} else {
				assert(valuetype == 1);

				if (src[i] != '{') {
					fprintf(stderr, "Error: Missing '{' after property name (line %d)\n", line);
					return false;
				}
				i++;

				float temp[3];
				for (int j = 0; j < 3; j++) {

					while (i < len && is_space(src[i])) {
						if (src[i] == '\n') line++;
						i++;
					}

					int sign = 1;
					if (src[i] == '-') {
						sign = -1;
						i++;
						if (i == len || !is_digit(src[i])) {
							fprintf(stderr, "Error: Missing number after minus sign (line %d)\n", line);
							return false;
						}
					} else if (!is_digit(src[i])) {
						fprintf(stderr, "Error: Missing number %d in vector value (line %d)\n", j, line);
						return false;
					}
					temp[j] = 0;
					do {
						int d = src[i] - '0';
						temp[j] = temp[j] * 10 + d;
						i++;
					} while (i < len && is_digit(src[i]));
					if (i < len && src[i] == '.') {
						i++; // Skip the dot
						if (i == len || !is_digit(src[i])) {
							fprintf(stderr, "Error: Missing decimal part after dot (line %d)\n", line);
							return false;
						}
						float q = 1.0f / 10;
						do {
							int d = src[i] - '0';
							temp[j] += q * d;
							q /= 10;
							i++;
						} while (i < len && is_digit(src[i]));
					}
					temp[j] *= sign;
				}

				while (i < len && is_space(src[i])) {
					if (src[i] == '\n') line++;
					i++;
				}

				if (i == len || src[i] != '}') {
					fprintf(stderr, "Error: Missing '}' after property value (line %d)\n", line);
					return false;
				}
				i++;

				value1.x = temp[0];
				value1.y = temp[1];
				value1.z = temp[2];
			}

			switch (prop) {

				case PROP_ALBEDO:
				if (value1.x < 0 || value1.x > 1 ||
					value1.y < 0 || value1.y > 1 ||
					value1.z < 0 || value1.z > 1) {
					fprintf(stderr, "Error: albedo values must be between 0 and 1 (line %d)\n", line);
					return false;
				}
				object.material.albedo = value1;
				break;

				case PROP_ROUGHNESS:
				if (value0 < 0 || value0 > 1) {
					fprintf(stderr, "Error: Roughness must be between 0 and 1 (line %d)\n", line);
					return false;
				}
				object.material.roughness = value0;
				break;

				case PROP_REFLECTANCE:
				if (value0 < 0 || value0 > 1) {
					fprintf(stderr, "Error: Reflectance must be between 0 and 1 (line %d)\n", line);
					return false;
				}
				object.material.reflectance = value0;
				break;

				case PROP_METALLIC:
				if (value0 < 0 || value0 > 1) {
					fprintf(stderr, "Error: Metallic must be between 0 and 1 (line %d)\n", line);
					return false;
				}
				object.material.metallic = value0;
				break;

				case PROP_EMISSION_POWER:
				object.material.emission_power = value0;
				break;

				case PROP_EMISSION_COLOR:
				if (value1.x < 0 || value1.x > 1 ||
					value1.y < 0 || value1.y > 1 ||
					value1.z < 0 || value1.z > 1) {
					fprintf(stderr, "Error: Emission color values must be between 0 and 1 (line %d)\n", line);
					return false;
				}
				object.material.emission_color = value1;
				break;

				case PROP_RADIUS:
				object.sphere.radius = value0;
				break;

				case PROP_CENTER:
				object.sphere.center = value1;
				break;

				case PROP_ORIGIN:
				object.cube.origin = value1;
				break;

				case PROP_SIZE:
				if (value1.x < 0 || value1.y < 0 || value1.z < 0) {
					fprintf(stderr, "Error: Size values must be positive (line %d)\n", line);
					return false;
				}
				object.cube.size = value1;
				break;
			}
		}

		if (scene->num_objects == MAX_OBJECTS)
			fprintf(stderr, "Warning: Ignoring object because the scene is too big (line %d)\n", line);
		else
			scene->objects[scene->num_objects++] = object;
	}

	return true;
}

bool parse_scene_file(char *file, Scene *scene)
{
	size_t len;
	char  *src = load_file(file, &len);
	if (src == NULL) {
		fprintf(stderr, "Error: Couldn't open scene file\n");
		return false;
	}

	bool ok = parse_scene_string(src, len, scene);

	free(src);
	return ok;
}
