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
	= "stroke-width=\"0.1px\" fill=\"none\"";

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
	const char * color = "#FF0000";

	printf("<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%s\" %s/>\n",
		x1, y1,
		x2, y2,
		color,
		stroke_string
	);
}


static void
svg_circle(
	const double x,
	const double y,
	const double rad,
	const char * const color
)
{
	printf("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"%s\" %s/>\n",
		x,
		y,
		rad,
		color,
		stroke_string
	);
}


int
main(void)
{
	stl_3d_t * const stl = stl_3d_parse(STDIN_FILENO);
	if (!stl)
		return EXIT_FAILURE;
	const double inset_distance = 5;
	const double hole_radius = 1.5;

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
		const int vertex_count = stl_trace_face(
			stl,
			f,
			vertex_list,
			face_used,
			0
		);

		fprintf(stderr, "%d: %d vertices\n", i, vertex_count);

		// generate a refernce frame based on this face
		refframe_t ref;
		refframe_init(&ref,
			f->vertex[0]->p,
			f->vertex[1]->p,
			f->vertex[2]->p
		);
		printf("<!-- face %d --><g>\n", i);

		// generate the polygon outline (should be one path?)
		for (int j = 0 ; j < vertex_count ; j++)
			svg_line(
				&ref,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p
			);

		// generate the inset mounting holes
		for (int j = 0 ; j < vertex_count ; j++)
		{
			double x, y;
			refframe_inset(&ref, inset_distance, &x, &y,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p,
				vertex_list[(j+2) % vertex_count]->p
			);
			svg_circle(x, y, hole_radius, "#00ff00");
		}

		printf("</g>\n");
	}

	printf("</g></svg>\n");

	return 0;
}
