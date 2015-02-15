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
	const double offset = 8;

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	const stl_vertex_t ** const vertex_list = calloc(sizeof(**vertex_list), stl->num_vertex);

	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];
		const v3_t origin = v->p;

		printf("// vertex %d\n"
			"translate([%f,%f,%f])\n"
			"//render() difference()\n"
			"{\n"
			"sphere(r=20);\n",
			i, origin.p[0], origin.p[1], origin.p[2]);
		
		int * const face_used = calloc(sizeof(*face_used), stl->num_face);

		for (int j = 0 ; j < v->num_face; j++)
		{
			// generate the polygon face for this vertex
			const stl_face_t * const f = v->face[j];
			if (face_used[f - stl->face])
				continue;
			const int start_vertex = v->face_num[j];
			const int vertex_count = stl_trace_face(
				stl,
				f,
				vertex_list,
				face_used,
				start_vertex
			);

			refframe_t ref;
			refframe_init(&ref,
				f->vertex[(start_vertex+0) % 3]->p,
				f->vertex[(start_vertex+1) % 3]->p,
				f->vertex[(start_vertex+2) % 3]->p
			);

			const stl_vertex_t * const v1 = vertex_list[0];
			const stl_vertex_t * const v2 = vertex_list[1];
			const v3_t d = v3_sub(v2->p, v1->p);
			const float len = v3_len(&v2->p, &v1->p);

			//const float b = acos(d.p[2] / len) * 180/M_PI;
			//const float c = d.p[0] == 0 ? sign(d.p[1]) * 90 : atan2(d.p[1], d.p[0]) * 180/M_PI;
//
			// use the transpose of the rotation matrix,
			// which will rotate from (x,y) to the correct
			// orientation relative to this connector node.

			printf("multmatrix(m=[[%f,%f,%f,0],[%f,%f,%f,0],[%f,%f,%f,0],[0,0,0,1]]) translate([0,0,%f]) linear_extrude(height=%f) polygon(points=[\n",
				ref.x.p[0],
				ref.y.p[0],
				ref.z.p[0],
				ref.x.p[1],
				ref.y.p[1],
				ref.z.p[1],
				ref.x.p[2],
				ref.y.p[2],
				ref.z.p[2],
				//a, b, c,
				-thickness/2,
				thickness
			);

			for(int k=0 ; k < vertex_count ; k++)
			{
				double x, y;
				v3_project(&ref, vertex_list[k]->p, &x, &y);
				printf("[%f,%f],", x, y);
			}
			printf("\n]);\n");

			// generate a polyhedron that spans
			// the width of this coplanar thingy
#if 0
			v3_t v0 = v3_sub(f->vertex[0]->p, origin);
			v3_t v1 = v3_sub(f->vertex[1]->p, origin);
			v3_t v2 = v3_sub(f->vertex[2]->p, origin);
			v3_t n;

			// compute normal of the face
			if (v->face_num[j] == 0)
				n = v3_cross(v1, v2);
			else
			if (v->face_num[j] == 1)
				n = v3_cross(v2, v0);
			else
			if (v->face_num[j] == 2)
				n = v3_cross(v0, v1);

			n = v3_scale(n, (thickness+1)/v3_mag(n)/2);

			// slide the vectors towards the center
			v3_t v0mid = v3_scale(v3_mid(v0, v1, v2), offset);
			v3_t v1mid = v3_scale(v3_mid(v1, v0, v2), offset);
			v3_t v2mid = v3_scale(v3_mid(v2, v0, v1), offset);

			v0 = v3_add(v0, v0mid);
			v1 = v3_add(v1, v1mid);
			v2 = v3_add(v2, v2mid);

			// compute the 
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
				"], %s = ["
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
				v5.p[0], v5.p[1], v5.p[2],
#ifdef __linux__
				"triangles"
#else
				"faces"
#endif
			);
#endif
				
			//break; // only do one right now
		}

		free(face_used);

		printf("}\n");
		//if (i == 2) break; // only do one right now
	}

	return 0;
}
