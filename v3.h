/** \file
 * 3D point operations.
 */
#ifndef _papercraft_v3_h_
#define _papercraft_v3_h_

#include <math.h>

#define EPS 0.001

#ifndef M_PI
#define 	M_PI   3.1415926535897932384
#endif

static inline float
sign(
	const float x
)
{
	if (x < 0)
		return -1;
	if (x > 0)
		return +1;
	return 0;
}


static inline float
min(
	const float a,
	const float b
)
{
	return a < b ? a : b;
}

static inline float
max(
	const float a,
	const float b
)
{
	return a > b ? a : b;
}

typedef struct
{
	float p[3];
} v3_t;


static inline int
v3_eq(
	const v3_t * v1,
	const v3_t * v2
)
{
	float dx = v1->p[0] - v2->p[0];
	float dy = v1->p[1] - v2->p[1];
	float dz = v1->p[2] - v2->p[2];

	if (-EPS < dx && dx < EPS
	&&  -EPS < dy && dy < EPS
	&&  -EPS < dz && dz < EPS)
		return 1;

	return 0;
}


static inline double
v3_len(
	const v3_t * const v0,
	const v3_t * const v1
)
{
	float dx = v0->p[0] - v1->p[0];
	float dy = v0->p[1] - v1->p[1];
	float dz = v0->p[2] - v1->p[2];

	return sqrt(dx*dx + dy*dy + dz*dz);
}

static inline double
v3_mag2(
	const v3_t v0
)
{
	float dx = v0.p[0];
	float dy = v0.p[1];
	float dz = v0.p[2];

	return dx*dx + dy*dy + dz*dz;
}

static inline double
v3_mag(
	const v3_t v0
)
{
	return sqrt(v3_mag2(v0));
}


static inline v3_t
v3_add(
	v3_t a,
	v3_t b
)
{
	v3_t c = { .p = {
		a.p[0] + b.p[0],
		a.p[1] + b.p[1],
		a.p[2] + b.p[2],
	} };
	return c;
}


static inline v3_t
v3_sub(
	v3_t a,
	v3_t b
)
{
	v3_t c = { .p = {
		a.p[0] - b.p[0],
		a.p[1] - b.p[1],
		a.p[2] - b.p[2],
	} };
	return c;
}

static inline v3_t
v3_scale(
	v3_t a,
	float s
)
{
	v3_t c = { .p = {
		a.p[0]*s,
		a.p[1]*s,
		a.p[2]*s,
	} };
	return c;
}


static inline
v3_t
v3_norm(
	const v3_t v
)
{
	return v3_scale(v, 1/v3_mag(v));
}


static inline
v3_t
v3_mid(
	const v3_t v0,
	const v3_t v1,
	const v3_t v2
)
{
	return v3_norm(
		v3_add(
			v3_sub(v1, v0),
			v3_sub(v2, v0)
		)
	);
}



static inline float
v3_dot(
	v3_t a,
	v3_t b
)
{
	return a.p[0]*b.p[0] + a.p[1]*b.p[1] + a.p[2]*b.p[2];
}


static inline v3_t
v3_cross(
	v3_t u,
	v3_t v
)
{
	float u1 = u.p[0];
	float u2 = u.p[1];
	float u3 = u.p[2];

	float v1 = v.p[0];
	float v2 = v.p[1];
	float v3 = v.p[2];

	v3_t c = { .p = {
		u2*v3 - u3*v2,
		u3*v1 - u1*v3,
		u1*v2 - u2*v1,
	}};

	return c;
}


// Compute the length of a line in screen space, ignoring Z
static inline float
v3_dist_2d(
	const v3_t * p0,
	const v3_t * p1
)
{
	const float dx = p1->p[0] - p0->p[0];
	const float dy = p1->p[1] - p0->p[1];

	return sqrt(dx*dx + dy*dy);
}


#endif
