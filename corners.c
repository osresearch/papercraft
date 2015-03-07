/** \file
 * Generate an OpenSCAD with connectors for each face.
 *
 * Options are inside only (with face flush on outside)
 * or with a slot for the face (like a corner cap)
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


static v3_t avg_x, avg_y, avg_z;

static void
print_multmatrix(
	const refframe_t * const ref,
	const int transpose
)
{
	printf("multmatrix(m=["
		"[%f,%f,%f,0],"
		"[%f,%f,%f,0],"
		"[%f,%f,%f,0],"
		"[ 0, 0, 0,1]])\n",
		transpose ? ref->x.p[0] : ref->x.p[0],
		transpose ? ref->x.p[1] : ref->y.p[0],
		transpose ? ref->x.p[2] : ref->z.p[0],
		transpose ? ref->y.p[0] : ref->x.p[1],
		transpose ? ref->y.p[1] : ref->y.p[1],
		transpose ? ref->y.p[2] : ref->z.p[1],
		transpose ? ref->z.p[0] : ref->x.p[2],
		transpose ? ref->z.p[1] : ref->y.p[2],
		transpose ? ref->z.p[2] : ref->z.p[2]
	);
}

static void
make_faces(
	const stl_3d_t * const stl,
	const stl_vertex_t * const v,
	const double thickness,
	const double translate,
	const double inset_dist,
	const double hole_dist,
	const double hole_rad,
	const double hole_height
)
{
	int * const face_used = calloc(sizeof(*face_used), stl->num_face);

	// generate all of the coplanar polygons at this vertex
	const stl_vertex_t ** const vertex_list = calloc(sizeof(**vertex_list), stl->num_vertex);

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

		avg_x = v3_add(avg_x, ref.x);
		avg_y = v3_add(avg_y, ref.y);
		avg_z = v3_add(avg_z, ref.z);

		// use the transpose of the rotation matrix,
		// which will rotate from (x,y) to the correct
		// orientation relative to this connector node.
		print_multmatrix(&ref, 0);
		printf("{\n");

		// generate the polygon plane
		if (thickness != 0)
		{
		printf("translate([0,0,%f]) linear_extrude(height=%f) polygon(points=[\n",
			translate,
			thickness
		);

		for(int k=0 ; k < vertex_count ; k++)
		{
			double x, y;
			refframe_inset(&ref, inset_dist, &x, &y,
				vertex_list[(k+0) % vertex_count]->p,
				vertex_list[(k+1) % vertex_count]->p,
				vertex_list[(k+2) % vertex_count]->p
			);
			printf("[%f,%f],", x, y);
		}
		printf("\n]);\n");
		}

		// generate the mounting holes/pins
		if (hole_rad != 0)
		{
		for(int k=0 ; k < vertex_count ; k++)
		{
			double x, y;
			refframe_inset(&ref, inset_dist+hole_dist, &x, &y,
				vertex_list[(k+0) % vertex_count]->p,
				vertex_list[(k+1) % vertex_count]->p,
				vertex_list[(k+2) % vertex_count]->p
			);
			printf("translate([%f,%f,0]) cylinder(r=%f,h=%f, $fs=1);\n",
				x, y,
				hole_rad,
				hole_height
			);
		}
		}

		printf("}\n");
	}

	free(face_used);
	free(vertex_list);
}


int
main(void)
{
	stl_3d_t * const stl = stl_3d_parse(STDIN_FILENO);
	if (!stl)
		return EXIT_FAILURE;
	const double thickness = 3;
	const double inset_dist = 2;
	const double hole_dist = 5;
	const double hole_rad = 3.3/2;

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits

	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];
		const v3_t origin = v->p;

		printf("// vertex %d\n"
			"//translate([%f,%f,%f])\n"
			"{\n"
			"render() difference()\n"
			"{\n",
			i, origin.p[0], origin.p[1], origin.p[2]);
		avg_x.p[0] = avg_x.p[1] = avg_x.p[2] = 0;
		avg_y.p[0] = avg_y.p[1] = avg_y.p[2] = 0;
		avg_z.p[0] = avg_z.p[1] = avg_z.p[2] = 0;
		
		//printf("render() intersection() {\n");
		printf("union() {\n");
		make_faces(stl, v, thickness, -thickness, 0, 0, 0, 0);

		printf("}\n");
		//printf("}\n");

		printf("union() {\n");
		// slice away the outer bits
		make_faces(stl, v, thickness, 0, -thickness, 0, 0, 0);
		//make_faces(stl, v, thickness, inset_dist, hole_dist, hole_rad);
		printf("}\n");


		printf("} // difference\n");

		// add back in the mounting pegs
		make_faces(stl, v, 0, 0, 0, hole_dist, hole_rad, thickness);
		printf("} // union\n");

		break;
		//if (i == 0) break; // only do one right now
	}

	refframe_t avg;
	refframe_init(&avg, avg_x, avg_y, avg_z);
	printf("%%");
	print_multmatrix(&avg, 1);
	//printf("translate([0,0,20]) sphere(r=2);\n");
	printf("cube([50,50,50], center=true);\n");
	return 0;
}
