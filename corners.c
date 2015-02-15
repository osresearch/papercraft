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
	const double thickness = 6;

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];
		const v3_t origin = v->p;

		printf("// vertex %d\n"
			"//translate([%f,%f,%f])\n"
			"//render() difference()\n"
			"{\n"
			"//sphere(r=20);\n",
			i, origin.p[0], origin.p[1], origin.p[2]);
		

		for (int j = 0 ; j < v->num_face; j++)
		{
			const stl_face_t * const f = v->face[j];
			v3_t v0 = v3_sub(f->vertex[0]->p, origin);
			v3_t v1 = v3_sub(f->vertex[1]->p, origin);
			v3_t v2 = v3_sub(f->vertex[2]->p, origin);
			v3_t n;

			if (v->face_num[j] == 0)
				n = v3_cross(v1, v2);
			else
			if (v->face_num[j] == 1)
				n = v3_cross(v2, v0);
			else
			if (v->face_num[j] == 2)
				n = v3_cross(v0, v1);

			n = v3_scale(n, (thickness+1)/v3_mag(n)/2);
			v3_t n2 = v3_scale(n, 1/v3_mag(n));

			v3_t v3 = v3_add(v0, n);
			v3_t v4 = v3_add(v1, n);
			v3_t v5 = v3_add(v2, n);
			v0 = v3_sub(v0, n);
			v1 = v3_sub(v1, n);
			v2 = v3_sub(v2, n);

			printf("polyhedron(\n"
				"points=[\n"
				"[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],\n"
				"[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],\n"
				"], faces = ["
				"  [0,1,2], [3,5,4],"
				"  [0,2,3], [2,5,3],"
				"  [0,3,4], [0,4,1],"
				"  [1,4,5], [1,5,2],"
				"]);\n",
				v0.p[0], v0.p[1], v0.p[2],
				v1.p[0], v1.p[1], v1.p[2],
				v2.p[0], v2.p[1], v2.p[2],
				v3.p[0], v3.p[1], v3.p[2],
				v4.p[0], v4.p[1], v4.p[2],
				v5.p[0], v5.p[1], v5.p[2]
			);
				
			//break; // only do one right now
		}

		printf("}\n");
		if (i == 0) break; // only do one right now
	}

	return 0;
}
