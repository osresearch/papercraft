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
	v3_t up
)
{
	camera_t * c = calloc(1, sizeof(*c));
	if (!c)
		return NULL;

	camera_setup(c, eye, lookat, up);
	return c;
}


void
camera_setup(
	camera_t * const c,
	v3_t eye,
	v3_t lookat,
	v3_t up
)
{
	// compute the basis for the camera
	// negative look direction from eye to destination
	v3_t w = v3_norm(v3_sub(eye, lookat));

	// compute the side axis
	v3_t u = v3_norm(v3_cross(up, w));

	// and the "up" normal
	v3_t v = v3_norm(v3_cross(w, u));

	float cam[4][4] = {};

	cam[0][0] = u.p[0];
	cam[1][0] = u.p[1];
	cam[2][0] = u.p[2];
	cam[3][0] = 0;

	cam[0][1] = v.p[0];
	cam[1][1] = v.p[1];
	cam[2][1] = v.p[2];
	cam[3][1] = 0;

	cam[0][2] = w.p[0];
	cam[1][2] = w.p[1];
	cam[2][2] = w.p[2];
	cam[3][2] = 0;

	// compute u dot c, v dot c, w, dot c
	cam[0][3] = -v3_dot(lookat, u);
	cam[1][3] = -v3_dot(lookat, v);
	cam[2][3] = -v3_dot(lookat, w);
	cam[3][3] = 1;

	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
			fprintf(stderr, " %+5.3f", cam[i][j]);
		fprintf(stderr, "\n");
	}

	// now compute the perspective projection matrix
	float fov = 60;
	float s = 1.0 / tan(fov * M_PI / 180 / 2);
	c->near = 1.0;
	c->far = 200;
	float f1 = - c->far * (c->far - c->near);
	float f2 = - c->far * c->near / (c->far - c->near);

	float pers[4][4] = {
		{ s, 0, 0, 0 },
		{ 0, s, 0, 0 },
		{ 0, 0, f1, -1 },
		{ 0, 0, f2, 0 },
	};

	// and apply it to the camera matrix to generate transform
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
		{
			float d = 0;
			for(int k = 0 ; k < 4 ; k++)
				d += cam[i][k] * pers[k][j];
			c->r[i][j] = d;
		}
	}
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
			p[i] += c->r[j][i] * v[j];

	if(0) fprintf(stderr, "%.2f,%.2f,%.2f -> %.2f,%.2f,%.2f,%.2f\n",
		v[0], v[1], v[2],
		p[0], p[1], p[2], p[3]
	);

/*
	for (int i = 0 ; i < 3 ; i++)
		p[i] /= p[3];
*/
	p[2] /= p[3];


	// Transform to screen coordinate frame,
	// and return it to the caller
	v_out->p[0] = p[0] * 100;
	v_out->p[1] = p[1] * 100;
	v_out->p[2] = p[2] * 100;

	// what if p->p[4] == 0?
	// pz < 0 == The point is behind us; do not display?
	if (p[2] < c->near || p[2] > c->far)
		return 0;

	return 1;
}
