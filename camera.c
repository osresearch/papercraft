/** \file
 * 3D camera equation.
 *
 * Given a camera matrix and a XYZ point, returns the 2D coordinate
 * of the point.
 */
#include <stdio.h>
#include <stdlib.h>

#include "camera.h"

struct _camera_t
{
	float	near;
	float	far;
	m44_t	r;
};


camera_t *
camera_new(
	v3_t eye,
	v3_t lookat,
	v3_t up,
	float fov
)
{
	camera_t * c = calloc(1, sizeof(*c));
	if (!c)
		return NULL;

	camera_setup(c, eye, lookat, up, fov);
	return c;
}


void
camera_setup(
	camera_t * const c,
	v3_t eye,
	v3_t lookat,
	v3_t up,
	float fov
)
{
	// compute the basis for the camera
	// negative look direction from eye to destination
	v3_t w = v3_norm(v3_sub(eye, lookat));

	// compute the side axis
	v3_t u = v3_norm(v3_cross(up, w));

	// and the "up" normal
	v3_t v = v3_norm(v3_cross(w, u));

	m44_t cam = {{
#if 0
		{ u.p[0], v.p[0], w.p[0], 0 },
		{ u.p[1], v.p[1], w.p[1], 0 },
		{ u.p[2], v.p[2], w.p[2], 0 },
		{ -v3_dot(u,eye), -v3_dot(v,eye), -v3_dot(w,eye), 1 },
#else
		{ u.p[0], u.p[1], u.p[2], -v3_dot(u,eye) },
		{ v.p[0], v.p[1], v.p[2], -v3_dot(v,eye) },
		{ w.p[0], w.p[1], w.p[2], -v3_dot(w,eye) },
		{ 0,      0,      0,      1 },
#endif
	}};


	fprintf(stderr, "Camera:\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", cam.m[i][j]);
		fprintf(stderr, "\n");
	}

	// now compute the perspective projection matrix
if(0) {
	float s = 1.0 / tan(fov * M_PI / 180 / 2);
	c->near = 1;
	c->far = 2;
	float f1 = - c->far / (c->far - c->near);
	float f2 = - c->far * c->near / (c->far - c->near);

	m44_t pers = {{
		{ s, 0, 0, 0 },
		{ 0, s, 0, 0 },
		{ 0, 0, f2, -1 },
		{ 0, 0, f1, 0 },
	}};

	fprintf(stderr, "Perspective:\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", pers.m[i][j]);
		fprintf(stderr, "\n");
	}
	// and apply it to the camera matrix to generate transform
	m44_mult(&c->r, &cam, &pers);
} else {
	// no perspective
	m44_t pers = {{
		{ 1, 0, 0, 0 },
		{ 0, 1, 0, 0 },
		{ 0, 0, 1, 0 },
		{ 0, 0, 0, 1 },
	}};
	// and apply it to the camera matrix to generate transform
	m44_mult(&c->r, &cam, &pers);
}


	fprintf(stderr, "Cam*Pers\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", c->r.m[i][j]);
		fprintf(stderr, "\n");
	}

}


/** Transform a XYZ point into a screen point.
 *
 * Returns 0 if this is behind us.  Perhaps it should do a z buffer?
 * https://en.wikipedia.org/wiki/3D_projection
 */
int
camera_project(
	const camera_t * const c,
	const v3_t * const v_in,
	v3_t * const v_out
)
{
	v4_t v = {{ v_in->p[0], v_in->p[1], v_in->p[2], 1 }};
	v4_t p = m44_multv(&c->r, &v);

	p.p[2] *= -1;

	// what if p->p[4] == 0?
	// pz < 0 == The point is behind us; do not display?
	//if (p[2] < c->near || p[2] > c->far)
	if (p.p[2] <= 0)
		return 0;

	// shrink by the distance
	p.p[0] *= 1000 / p.p[2];
	p.p[1] *= 1000 / p.p[2];
	//p[2] /= 1000;

	// Transform to screen coordinate frame,
	// and return it to the caller
	v_out->p[0] = p.p[0];
	v_out->p[1] = p.p[1];
	v_out->p[2] = p.p[2];

	v3_print(*v_out);

	return 1;
}
