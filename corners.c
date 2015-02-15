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
		const v3_t origin = v->p;

		printf("// vertex %d\n"
			"translate([%f,%f,%f]) {\n"
			"sphere(r=20);\n",
			i, origin.p[0], origin.p[1], origin.p[2]);
		

		for (int j = 0 ; j < v->num_face; j++)
		{
			const stl_face_t * const f = v->face[j];
			v3_t v1 = v3_sub(f->vertex[0]->p, origin);
			v3_t v2 = v3_sub(f->vertex[1]->p, origin);
			v3_t v3 = v3_sub(f->vertex[2]->p, origin);

			printf("%%polyhedron(\n"
				"points=[\n"
				"[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],\n"
				"], faces = [[0,1,2]]);\n",
				v1.p[0], v1.p[1], v1.p[2],
				v2.p[0], v2.p[1], v2.p[2],
				v3.p[0], v3.p[1], v3.p[2]
			);
				
/*
			fprintf(stderr, "\t%d: %d:%f,%d:%f,%d:%f\n",
				f - stl->face,
				f->face[0] ? f->face[0] - stl->face : -1, f->angle[0],
				f->face[1] ? f->face[1] - stl->face : -1, f->angle[1],
				f->face[2] ? f->face[2] - stl->face : -1, f->angle[2]
			);
*/
		}

		printf("}\n");
		break;
	}

	return 0;
}
