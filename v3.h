/** \file
 * 3D point operations.
 */
#ifndef _papercraft_v3_h_
#define _papercraft_v3_h_

#include <math.h>

#define EPS 0.5


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

#endif
