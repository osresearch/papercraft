#pragma once

#include "v3.h"

typedef struct _camera_t camera_t;

extern camera_t *
camera_new(
	float eye_z,
	float phi,
	float theta,
	float psi
);

extern void
camera_setup(
	camera_t * const c,
	float eye_z,
	float phi,
	float theta,
	float psi
);

/** Transform a XYZ point into a screen point.
 *
 * Returns 0 if this is behind us.  Perhaps it should do a z buffer?
 */
extern int
camera_project(
	const camera_t * const c,
	const v3_t * const v,
	v3_t * const v_out
);


