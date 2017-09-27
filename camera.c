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
	float	zoom;
	v3_t	eye;
	float	r[3][3];
};


camera_t *
camera_new(
	v3_t eye,
	float phi,
	float theta,
	float psi
)
{
	camera_t * c = calloc(1, sizeof(*c));
	if (!c)
		return NULL;

	c->zoom = 4096;
	camera_setup(c, eye, phi, theta, psi);
	return c;
}


void
camera_setup(
	camera_t * const c,
	v3_t eye,
	float phi,
	float theta,
	float psi
)
{
	const float sx = sin(phi);
	const float cx = cos(phi);
	const float sy = sin(theta);
	const float cy = cos(theta);
	const float sz = sin(psi);
	const float cz = cos(psi);

	c->r[0][0] =   cy * cz;
	c->r[0][1] = (-cy * sz) + (sx * sy * cz);
	c->r[0][2] = ( sx * sz) + (cx * sy * cz);

	c->r[1][0] =   cx * sz;
	c->r[1][1] = ( cx * cz) + (sx * sy * sz);
	c->r[1][2] = (-sx * cz) + (cx * sy * sz);

	c->r[2][0] = -sy;
	c->r[2][1] =  sx * cy;
	c->r[2][2] =  cx * cy;

	c->eye = eye;
 }


/** Transform a XYZ point into a screen point.
 *
 * Returns 0 if this is behind us.  Perhaps it should do a z buffer?
 */
int
camera_project(
	const camera_t * const c,
	const v3_t * const v,
	v3_t * const v_out
)
{
	v3_t p = c->eye;

	for (int i = 0 ; i < 3 ; i++)
		for (int j = 0 ; j < 3 ; j++)
			p.p[i] += c->r[i][j] * v->p[j];

	// what if p->p[2] == 0?

	// Transform to screen coordinate frame,
	// is this rotating?
	float px = p.p[1] / p.p[2];
	float py = p.p[0] / p.p[2];
	float pz = p.p[2];

	// return it to the caller
	v_out->p[0] = px * c->zoom;
	v_out->p[1] = py * c->zoom;
	v_out->p[2] = pz * c->zoom;

	// pz < 0 == The point is behind us; do not display?
	return pz <= 0 ? 0 : 1;
}
