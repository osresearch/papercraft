/** \file
 * Generate an OpenSCAD with connectors for each face.
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

static int debug = 0;
static int draw_labels = 0;


int
main(void)
{
	stl_3d_t * const stl = stl_3d_parse(STDIN_FILENO);
	if (!stl)
		return EXIT_FAILURE;

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];

		fprintf(stderr, "%d: %f,%f,%f\n",
			i, v->p.p[0], v->p.p[1], v->p.p[2]);

		for (int j = 0 ; j < v->num_face; j++)
		{
			const stl_face_t * const f = v->face[j];

			fprintf(stderr, "\t%d: %d:%f,%d:%f,%d:%f\n",
				f - stl->face,
				f->face[0] ? f->face[0] - stl->face : -1, f->angle[0],
				f->face[1] ? f->face[1] - stl->face : -1, f->angle[1],
				f->face[2] ? f->face[2] - stl->face : -1, f->angle[2]
			);
		}
	}
#if 0
		if (debug) fprintf(stderr, "---------- triangle %d (%d)\n", i, num_vertex);

		stl_vertex_t * vp[3] = {};

		for (int j = 0 ; j < 3 ; j++)
		{
			const v3_t * const p = &stl_faces[i].p[j];
			vp[j] = stl_vertex_find(vertices, &num_vertex, p);
		}

		// walk all of other triangles to figure out if
		// any of the triangles are coplanar and have shared
		// edges.

		uint8_t coplanar_mask = 0;
		for (int j = 0 ; j < num_triangles ; j++)
		{
			if (j == i)
				continue;

			if (debug)
			fprintf(stderr, "%d: check %d -> %d\n", num_vertex, i, j);

			coplanar_mask |= coplanar_check(
				&stl_faces[i], &stl_faces[j]);
		}

		if (debug)
			fprintf(stderr, "mask %d\n", coplanar_mask);

		// all three vertices are mapped; generate the
		// connections
		for (int j = 0 ; j < 3 ; j++)
		{
			stl_vertex_t * const v = vp[j];

			// if the edge from j to j+1 is not coplanar,
			// add it to the list
			if ((coplanar_mask & (1 << j)) == 0)
			{
				if (debug)
				fprintf(stderr, "%p: %d insert\n", v, j);
				stl_edge_insert(v, vp[(j+1) % 3]);
			}
		}
	}

	fprintf(stderr, "%d unique vertices\n", num_vertex);
	printf("thick=%f;\n"
		"module connector(len) {\n"
		"  render() difference() {\n"
		"    cylinder(r=thick/2+2, h=2*thick);\n"
		//"    translate([0,0,len/2+2]) cube([thick,thick,2*thick]);\n"
		"    translate([0,0,thick/2+2]) cylinder(r=thick/2, h=2*thick);\n"
		"  }\n"
		//"  %%translate([0,0,len*0.48/2]) cube([thick,thick,len*0.48], center=true);\n"
		"  %%translate([0,0,0]) cylinder(r=thick/2, h=len*0.48);\n"
		"}\n",
		thick
	);

	for (int i = 0 ; i < num_vertex ; i++)
	{
		stl_vertex_t * const v = vertices[i];
		printf("translate([%f,%f,%f]) {\n",
			v->p.p[0],
			v->p.p[1],
			v->p.p[2]
		);
		
		printf("sphere(r=%f); // %d %p\n", thick/2+2, i, v);

		for (int j = 0 ; j < v->num_edges ; j++)
		{
			stl_vertex_t * const v2 = v->edges[j];
			const v3_t d = v3_sub(v2->p, v->p);
			const float len = v3_len(&v2->p, &v->p);

			const float b = acos(d.p[2] / len) * 180/M_PI;
			const float c = d.p[0] == 0 ? sign(d.p[1]) * 90 : atan2(d.p[1], d.p[0]) * 180/M_PI;
//
			printf("rotate([0,%f,%f]) ", b, c);

			if (do_square)
				printf("connector(%f);\n", len);
			else
				printf(" cylinder(r=1, h=%f); // %p\n",
					len*.45,
					v2
				);
		}

		printf("}\n");
	}
#endif

	return 0;
}
