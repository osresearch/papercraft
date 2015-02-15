/** \file
 * Generate an svg file with the polygonal faces.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <err.h>
#include <assert.h>
#include "v3.h"
#include "stl_3d.h"

static const char * stroke_string
	= "stroke=\"#FF0000\" stroke-width=\"1px\" fill=\"none\"";

typedef struct
{
	v3_t origin;
	v3_t x;
	v3_t y;
	v3_t z;
} refframe_t;


static void
v3_project(
	const refframe_t * const ref,
	const v3_t p_in,
	double * const x_out,
	double * const y_out
)
{
	v3_t p = v3_sub(p_in, ref->origin);

	double x = ref->x.p[0]*p.p[0] + ref->x.p[1]*p.p[1] + ref->x.p[2]*p.p[2];
	double y = ref->y.p[0]*p.p[0] + ref->y.p[1]*p.p[1] + ref->y.p[2]*p.p[2];
	double z = ref->z.p[0]*p.p[0] + ref->z.p[1]*p.p[1] + ref->z.p[2]*p.p[2];

	fprintf(stderr, "%f,%f,%f\n", x, y, z);

	*x_out = x;
	*y_out = y;
}


static void
svg_line(
	const refframe_t * const ref,
	const v3_t p1,
	const v3_t p2
)
{
	// project p1 and p2 into the plane
	double x1, y1, x2, y2;

	v3_project(ref, p1, &x1, &y1);
	v3_project(ref, p2, &x2, &y2);

	printf("<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" %s/>\n",
		x1, y1,
		x2, y2,
		stroke_string
	);
}


void
draw_face(
	const stl_3d_t * const stl,
	const stl_face_t * const f,
	int * const face_used,
	const refframe_t * const ref
)
{
	if (face_used[f - stl->face]++)
		return;

	for (int i = 0 ; i < 3 ; i++)
	{
		if (!f->face[i] || f->angle[i] != 0)
		{
			// not coplanar (or missing?)
			svg_line(ref, f->vertex[i]->p, f->vertex[(i+1)%3]->p);
		} else {
			// coplanar, recurse
			draw_face(stl, f->face[i], face_used, ref);
		}
	}
}


int
main(void)
{
	stl_3d_t * const stl = stl_3d_parse(STDIN_FILENO);
	if (!stl)
		return EXIT_FAILURE;
	const double inset = 8;
	const double hole = 3;

	int * const face_used = calloc(sizeof(*face_used), stl->num_face);

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	printf("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");
	printf("<g transform=\"scale(3.543307)\"><!-- scale to mm -->\n");


	for(int i = 0 ; i < stl->num_face ; i++)
	{
		if (face_used[i])
			continue;

		const stl_face_t * const f = &stl->face[i];

		// generate a refernce frame based on this face
		refframe_t ref = {
			.origin = f->vertex[0]->p,
		};

		const v3_t dx = v3_norm(v3_sub(f->vertex[1]->p, ref.origin));
		const v3_t dy = v3_norm(v3_sub(f->vertex[2]->p, ref.origin));

		ref.x = dx;
		ref.z = v3_norm(v3_cross(dx, dy));
		ref.y = v3_norm(v3_cross(ref.x, ref.z));

		printf("<!-- face %d --><g>\n", i);
		draw_face(stl, f, face_used, &ref);
		printf("</g>\n");
	}

	printf("</g></svg>\n");

	return 0;
}
