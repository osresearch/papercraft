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


/** Starting at a point, trace the coplanar polygon and return a
 * list of vertices.
 */
int
trace_face(
	const stl_3d_t * const stl,
	const stl_face_t * const f_start,
	const stl_vertex_t ** vertex_list,
	int * const face_used
)
{
	const stl_face_t * f = f_start;
	int i = 0;
	int vertex_count = 0;

	// preload the list with the starting point
	vertex_list[vertex_count++] = f->vertex[i];

	do {
		const stl_vertex_t * const v1 = f->vertex[(i+0) % 3];
		const stl_vertex_t * const v2 = f->vertex[(i+1) % 3];
		const stl_face_t * const f_next = f->face[i];

		fprintf(stderr, "%p %d: %f,%f,%f\n", f, i, v1->p.p[0], v1->p.p[1], v1->p.p[2]);
		face_used[f - stl->face] = 1;

		if (!f_next || f->angle[i] != 0)
		{
			// not coplanar or no connection.
			// add the NEXT vertex on this face and continue
			vertex_list[vertex_count++] = v2;
			i = (i+1) % 3;
			continue;
		}

		// coplanar; figure out which vertex on the next
		// face we start at
		int i_next = -1;
		for (int j = 0 ; j < 3 ; j++)
		{
			if (f_next->vertex[j] != v1)
				continue;
			i_next = j;
			break;
		}

		if (i_next == -1)
			abort();
		
		// move to the new face
		f = f_next;
		i = i_next;

		// keep going until we reach our starting face
		// at vertex 0.
	} while (f != f_start || i != 0);

	return vertex_count;
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

	const stl_vertex_t ** const vertex_list = calloc(sizeof(*vertex_list), stl->num_vertex);

	for(int i = 0 ; i < stl->num_face ; i++)
	{
		if (face_used[i])
			continue;

		const stl_face_t * const f = &stl->face[i];
		const int vertex_count = trace_face(stl, f, vertex_list, face_used);

		fprintf(stderr, "%d: %d vertices\n", i, vertex_count);

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

		for (int j = 0 ; j < vertex_count-1 ; j++)
			svg_line(&ref, vertex_list[j]->p, vertex_list[j+1]->p);
		printf("</g>\n");
	}

	printf("</g></svg>\n");

	return 0;
}
