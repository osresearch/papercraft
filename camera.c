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
	float	r[4][4];
};


camera_t *
camera_new(
	v3_t eye,
	v3_t lookat,
	v3_t up,
	float fov,
	float scale
)
{
	camera_t * c = calloc(1, sizeof(*c));
	if (!c)
		return NULL;

	camera_setup(c, eye, lookat, up, fov, scale);
	return c;
}


void
camera_setup(
	camera_t * const c,
	v3_t eye,
	v3_t lookat,
	v3_t up,
	float fov,
	float scale
)
{
	// compute the basis for the camera
	// negative look direction from eye to destination
	v3_t w = v3_norm(v3_sub(eye, lookat));

	// compute the side axis
	v3_t u = v3_norm(v3_cross(up, w));

	// and the "up" normal
	v3_t v = v3_norm(v3_cross(w, u));

	float cam[4][4] = {
#if 0
		{ u.p[0], v.p[0], w.p[0], 0 },
		{ u.p[1], v.p[1], w.p[1], 0 },
		{ u.p[2], v.p[2], w.p[2], 0 },
		{ -v3_dot(u,eye), -v3_dot(v,eye), -v3_dot(w,eye), 1 },
#else
		{ u.p[0], u.p[1], u.p[2], -v3_dot(u,lookat) },
		{ v.p[0], v.p[1], v.p[2], -v3_dot(v,lookat) },
		{ w.p[0], w.p[1], w.p[2], -v3_dot(w,lookat) },
		{ 0,      0,      0,      1 },
#endif
	};


	fprintf(stderr, "Camera:\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", cam[i][j]);
		fprintf(stderr, "\n");
	}

	// now compute the perspective projection matrix
	float s = 1.0 / tan(fov * M_PI / 180 / 2);
	c->near = 1.0;
	c->far = 20;
	float f1 = - c->far / (c->far - c->near);
	float f2 = - c->far * c->near / (c->far - c->near);

	float pers[4][4] = {
		{ s, 0, 0, 0 },
		{ 0, s, 0, 0 },
		{ 0, 0, f1, -1 },
		{ 0, 0, f2, 0 },
	};

	fprintf(stderr, "Perspective:\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", pers[i][j]);
		fprintf(stderr, "\n");
	}

	// and apply it to the camera matrix to generate transform
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
		{
			float d = 0;
			for(int k = 0 ; k < 4 ; k++)
				d += pers[i][k] * cam[k][j];
			c->r[i][j] = d;
		}
	}


	fprintf(stderr, "Cam*Pers\n");
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", c->r[i][j]);
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
	float v[4] = { v_in->p[0], v_in->p[1], v_in->p[2], 1 };
	float p[4] = { 0, 0, 0, 0};

	for (int i = 0 ; i < 4 ; i++)
		for (int j = 0 ; j < 4 ; j++)
			p[i] += c->r[i][j] * v[j];

	if(1) fprintf(stderr, "%.2f %.2f %.2f -> %.2f %.2f %.2f %.2f\n",
		v[0], v[1], v[2],
		p[0], p[1], p[2], p[3]
	);

/*
	for (int i = 0 ; i < 3 ; i++)
		p[i] /= p[3];
	p[0] *= -1;
	p[1] *= -1;
	p[2] /= p[3];
*/


	// Transform to screen coordinate frame,
	// and return it to the caller
	v_out->p[0] = p[0];
	v_out->p[1] = p[1];
	v_out->p[2] = p[2];

	// what if p->p[4] == 0?
	// pz < 0 == The point is behind us; do not display?
	if (p[2] < c->near || p[2] > c->far)
		return 0;

	return 1;
}
