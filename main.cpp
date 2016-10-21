#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <array>


enum {
#ifdef NICE
	WIDTH	= 1600,
	HEIGHT	= 1200,
	SAMPLES = 400,
#else
	WIDTH	= 800,
	HEIGHT	= 600,
	SAMPLES = 20,
#endif
};

const float E = 0.0005;

struct vec {
	float x, y, z;
	vec(float x=0, float y=0, float z=0) : x(x), y(y), z(z) {}
	vec operator +(vec const & o) { return vec(x+o.x, y+o.y, z+o.z); }
	vec operator -(vec const & o) { return vec(x-o.x, y-o.y, z-o.z); }
	vec operator *(vec const & o) { return vec(x*o.x, y*o.y, z*o.z); }
	vec operator *(float f) { return vec(x*f, y*f, z*f); }
	float length() { return sqrtf(dot(*this)); }
	float lengthsq() { return dot(*this); }
	vec normalized() { return *this * (1.0 / length()); }
	float dot(vec o) { return x*o.x + y*o.y + z*o.z; }
	vec cross(vec o) { return vec(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x); }
};

float sphere(vec p, float r=1) {
	return p.length() - r;
}

float box(vec p, vec s=vec(1, 1, 1), float r=0) {
	vec d = vec(fabs(p.x), fabs(p.y), fabs(p.z)) - s;
	return vec(std::max(d.x, 0.0f), std::max(d.y, 0.0f), std::max(d.z, 0.0f)).length() - r;
}

float ring(vec p, float r1=1, float r2=0.2) {
	return vec(vec(p.x, p.y).length() - r1, p.z).length() - r2;
}

float cylinder(vec p, float r=0.2, float h=1) {
	return std::max(vec(p.x, p.z).length() - r, fabsf(p.y) - h);
}

enum Material {
	MAT_FLOOR,
	MAT_FLOOR2,
	MAT_SPHERE,
	MAT_BOX,
	MAT_BOX2,
	MAT_BOX3,
};

vec colors[] = {
	{ 0.9, 0.8, 0.8 },
	{ 0.3, 0.3, 0.3 },
	{ 0.9, 0.8, 0.4 },
	{ 0.9, 0.4, 0.4 },
	{ 0.5, 0.8, 1 },
	{ 0.7, 0.8, 0.5 },
};

// material of current object
Material material;

float combine(float a, float b, Material mat) {
	if (b < a) material = mat;
	return std::min(a, b);
}

float combine_round(float a, float b, float r=0.1) {
	if (std::max(a, b) < r) return r - vec(r - a, b - r).length();
	return std::min(a, b);
}

float smin(float a, float b, float k=0.2) {
	float h = std::max(0.0f, std::min(1.0f, 0.5f + 0.5f * (b-a) / k));
	return std::min(h, std::min(a, b)) - h * (1 - h) * k;
}

vec rotate_y(vec p, float ang) {
	float s = sinf(ang);
	float c = cosf(ang);
	return vec(p.x * c + p.z * s, p.y, p.z * c - p.x * s);
}

float frand() { return rand() * 2.0 / float(RAND_MAX) - 1; }

float distance(vec p) {
	material = MAT_FLOOR;
	float d = p.y + 1;

	vec ps;

	ps = p + vec(1.5, 0, 0.2);
	d = combine(d,
		combine_round(
			std::max(sphere(ps, 1), ps.y),
			sphere(ps, 0.8), 0.2),
		MAT_SPHERE);

	ps = rotate_y(p, 0.1) - vec(1, -0.4, 0);
	d = combine(d,
		-combine_round(
			-box(ps, vec(0.6, 0.2, 0.6), 0.4),
			sphere(ps - vec(0, 0.7, 0), 0.7),
			0.1),
		MAT_BOX);

	vec pb = rotate_y(p, 0.3) - vec(0, 0.5, 2);
	pb.x = fabsf(pb.x);
	vec pc = pb - vec(0.5, 0, -0.25);
	d = combine(d,
		std::min(
			box(pb, vec(1, 1.5, 0.8)),
			box(pc, vec(0.4, 1.4, 0.6), 0.05)),
		MAT_BOX2);

	float s = std::max(0.0f, std::min(1.0f, p.y));
	float z = std::max(0.0f, std::min<float>(3.14159 * 6, p.y * 20));
	s = s * s * (3 - 2 * s);
	d = combine(d,
		box(rotate_y(p - vec(2.4, 0.3, 1.7), -0.5),
			vec(0.9 - s * 0.4, 1.1, 0.5 - s * 0.2), 0.2 + -cosf(z) * 0.03),
		MAT_BOX3);

	return d;
}

vec normal(vec p) {
	float d = distance(p);
	return vec(
		d - distance(p - vec(E, 0, 0)),
		d - distance(p - vec(0, E, 0)),
		d - distance(p - vec(0, 0, E))).normalized();
}

vec random_point_in_shpere() {
	vec v;
	do v = vec(frand(), frand(), frand());
	while (v.lengthsq() < 1);
	return v;
}

vec find_color(vec sp) {
	vec px = (sp - vec(WIDTH, HEIGHT) * 0.5) * (2.0 / WIDTH);
	vec dx = (px * 0.85 - vec(0, 0, -1)).normalized();

	float const MAX_DIST = 50;
	vec col_sum;
	for (int fn = 0; fn < SAMPLES; ++fn) {
		vec pos = vec(1.6, 2.4, -4) + random_point_in_shpere() * frand() * 0.04;
		vec look_at = vec(0.55, 0.3, 0);

		vec dir = (look_at - pos).normalized();
		vec left = dir.cross(vec(0, 1, 0));
		vec up = left.cross(dir);
		dir = dir * dx.z + left * dx.x + up * dx.y;

		float t = 0, t2 = 0;
		float factor = 0.9;
		int i;
		const int N = 240;
		vec p;
		for (i = 0; i < N && t <= MAX_DIST; ++i) {
			p = pos + dir * t;
			float d = distance(p);
			if (d < 0) {
				factor *= 0.9;
				t = t2;
				continue;
			}
			if (d < E * 0.1) break;
			t2 = t;
			t += std::min(d * factor, 1.0f);
		}

		vec col_fog = vec(0.15, 0.15, 0.25);
		vec col = col_fog;
		if (i < N && t <= MAX_DIST) {

			// checker board
			if (material == MAT_FLOOR
			&& (fmod(p.x + 100, 2) < 1) ^ (fmod(p.z + 100, 2) < 1)) {
				material = MAT_FLOOR2;
			}

			vec col_diffuse = colors[material];
			vec col_ambient = colors[material] * 0.9;
			vec col_specular;

			vec dir_light = vec(0.5, -1, 0.5).normalized();
			vec norm = normal(p);
//			if (0)
			{
				float f = std::max(0.0f, -norm.dot(dir_light));
				col_diffuse = col_diffuse * f;
			}
//			if (0)
			{
				vec reflect = dir_light + (norm * norm.dot(dir_light)) * -2;
				float f = pow(reflect.dot(dir), 20);
				col_specular = vec(1, 1, 1) * f;
			}
//			if (0)
			{
				const float SHADOW_DIST = 10;
				const float SHADOW_N = 100;
				const float SMOOTH_SHADOW = 0.2;
				float f = SMOOTH_SHADOW;
				float t = E;
				for (int i = 0; i < SHADOW_N && t < SHADOW_DIST; ++i) {
					float d = distance(p - dir_light * t);
					t += std::max(E * 0.1f, d * 0.7f);
					f = std::max(0.0f, std::min(f, d / t));
				}
				f /= SMOOTH_SHADOW;
				col_diffuse = col_diffuse * f;
				col_specular = col_specular * f;
			}
//			if (0)
			{
				const float AC_DIST = 2;
				vec q = p + norm * E * 2;
				vec r = random_point_in_shpere().normalized();
				if (r.dot(norm) < 0) r = r * -1;
				float t = 0;
				while (t < AC_DIST) {
					float d = distance(q + r * t);
					t += std::min(d * 0.7f, 0.3f);
					if (d < E) break;
				}
				t = std::min(t, AC_DIST);
				col_ambient = col_ambient * powf(t / AC_DIST, 2.5);
			}

			float f = pow(0.88, t); // fog
			col = (col_diffuse + col_ambient + col_specular) * f + col_fog * (1 - f);
		}
		col_sum = col_sum + col;
	}
	col_sum = col_sum * (1.0 / SAMPLES);
	col_sum = col_sum * powf((px + vec(0, 0, 3.5)).normalized().dot(vec(0, 0, 1.01)), 30);
	return col_sum;
}

int main() {
	FILE* f = fopen("out.ppm", "w");
	fprintf(f, "P6 %d %d 255\n", WIDTH, HEIGHT);
	for (int y = 0; y < HEIGHT; ++y) {
		printf("\rprogress:%3d%%", y * 100 / HEIGHT);
		fflush(stdout);
		for (int x = 0; x < WIDTH; ++x) {
			vec c = find_color(vec(x, HEIGHT - y));
			fprintf(f, "%c%c%c",
				std::max(0, std::min<int>(255, c.x * 255)),
				std::max(0, std::min<int>(255, c.y * 255)),
				std::max(0, std::min<int>(255, c.z * 255)));
		}
	}
	fclose(f);
	return 0;
}
