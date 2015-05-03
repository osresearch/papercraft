/** \file
 * Generate an OpenSCAD with connectors for each face.
 *
 * This imports the original STL file and then slices the corners
 * off from it. 
 * Options are inside only (with face flush on outside)
 * or with a slot for the face (like a corner cap)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
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
print_normal(
	const v3_t * normal
)
{
	const float x = normal->p[0];
	const float y = normal->p[1];
	const float z = normal->p[2];
	const double length = sqrt(x*x+y*y+z*z);
	const double b = acos(z / length);
	const double c = x == 0 ? sign(y)*90 : atan2(y,x);

	printf("rotate([%f,%f,%f])\n", 0.0, b * 180 / M_PI, c * 180 / M_PI);
}


static void
find_normal(
	const stl_3d_t * const stl,
	const stl_vertex_t * const v,
	const float inset_dist,
	v3_t * const avg
)
{
	int * const face_used
		= calloc(sizeof(*face_used), stl->num_face);

	// generate all of the coplanar polygons at this vertex
	const stl_vertex_t ** const vertex_list
		= calloc(sizeof(**vertex_list), stl->num_vertex);

	for (int j = 0 ; j < v->num_face; j++)
	{
		// generate the polygon face for this vertex
		const stl_face_t * const f = v->face[j];
		if (face_used[f - stl->face])
			continue;

		//ref.origin.p[0] = 0;
		//ref.origin.p[1] = 0;
		//ref.origin.p[2] = 0;

		const int start_vertex = v->face_num[j];
		const int vertex_count = stl_trace_face(
			stl,
			f,
			vertex_list,
			face_used,
			start_vertex
		);

		// find this vertex in the vertex list
		// and compute the vector that subdivides the
		// two outbound edges
		for (int k = 0 ; k < vertex_count ; k++)
		{
			if (vertex_list[k] != v)
				continue;

			v3_t p1 = vertex_list[(k+vertex_count-1) % vertex_count]->p;
			v3_t p2 = vertex_list[k % vertex_count]->p;
			v3_t p3 = vertex_list[(k+1) % vertex_count]->p;

			refframe_t ref;
			refframe_init(
				&ref,
				p2,
				p3,
				p1
			);

			double x, y;
			refframe_inset(&ref, inset_dist, &x, &y, p1, p2, p3);

			v3_t hole = refframe_project(&ref, (v3_t){{x,y,0}});
			//hole = refframe_project(&ref, (v3_t){{10,0,0}});
			//hole.p[0] = 10*ref.x.p[0]; // + ref.origin.p[0];
			//hole.p[1] = 10*ref.x.p[1]; // + ref.origin.p[1];
			//hole.p[2] = 10*ref.x.p[2]; // + ref.origin.p[2];
			fprintf(stderr, "**** %p [%f,%f]=>%f,%f,%f\n", v, x, y,
				hole.p[0],
				hole.p[1],
				hole.p[2]
			);

			printf("color(\"green\") translate([%f,%f,%f]) sphere(r=1);\n",
				hole.p[0],
				hole.p[1],
				hole.p[2]
			);

#if 0
			hole = refframe_project(&ref, (v3_t){10,0,0});
			//v3_t hole = refframe_project(&ref, (v3_t){5,5,0});
			printf("translate([%f,%f,%f]) sphere(r=1);\n",
				hole.p[0],
				hole.p[1],
				hole.p[2]
			);
/*
			hole = refframe_project(&ref, (v3_t){0,10,0});
			//v3_t hole = refframe_project(&ref, (v3_t){5,5,0});
			printf("%%translate([%f,%f,%f]) sphere(r=1);\n",
				hole.p[0],
				hole.p[1],
				hole.p[2]
			);
*/
#endif

			*avg = v3_add(*avg, ref.z);
			//*avg = v3_add(*avg, ref.x);

			//avg_x = v3_add(avg_x, ref.x);
			//avg_y = v3_add(avg_y, ref.y);
			//*avg = v3_add(*avg, p);
		}
		//return;
	}

#if 0
		// use the transpose of the rotation matrix,
		// which will rotate from (x,y) to the correct
		// orientation relative to this connector node.
		print_multmatrix(&ref, 1);
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
				printf("translate([%f,%f,%f]) cylinder(r=%f,h=%f, $fs=1);\n",
					x, y, -hole_height/2,
					hole_rad,
					hole_height
				);
			}
		}

		printf("}\n");
	}
#endif

	free(face_used);
	free(vertex_list);
}


int
main(
	int argc,
	char ** argv
)
{
	if (argc <= 1)
	{
		fprintf(stderr, "Usage: corners file.stl > file-corners.scad\n");
		return -1;
	}

	const char * const stl_name = argv[1];
	int fd = open(stl_name, O_RDONLY);
	if (fd < 0)
	{
		perror(stl_name);
		return -1;
	}

	stl_3d_t * const stl = stl_3d_parse(fd);
	if (!stl)
		return EXIT_FAILURE;
	close(fd);

	printf("module model() {\n"
		"render() difference() {\n"
	 	"import(\"%s\");\n",
		stl_name
	);
	//printf("%%model();\n");

	const double thickness = 3;
	const double inset_dist = 5;
	const double hole_dist = 5;
	const double hole_rad = 1.25;

	int * const face_used
		= calloc(sizeof(*face_used), stl->num_face);
	const stl_vertex_t ** const vertex_list
		= calloc(sizeof(*vertex_list), stl->num_vertex);

	// for face, generate the set of coplanar points that go with it
	// and "drill" holes in the model for those corners.
	for (int i = 0 ; i < stl->num_face ; i++)
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

		refframe_t ref;
		refframe_init(
			&ref,
			f->vertex[0]->p,
			f->vertex[1]->p,
			f->vertex[2]->p
		);

		// replace the origin with the actual origin
		//ref.origin.p[0] = 0;
		//ref.origin.p[1] = 0;
		//ref.origin.p[2] = 0;

		printf("translate([%f,%f,%f])",
			f->vertex[0]->p.p[0],
			f->vertex[0]->p.p[1],
			f->vertex[0]->p.p[2]
		);
			
		print_multmatrix(&ref, 0);
		printf("{\n");

		// generate a bolt hole for each non-copolanar corner
		for (int j = 0 ; j < vertex_count ; j++)
		{
			double x, y;
			refframe_inset(
				&ref,
				inset_dist,
				&x,
				&y,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p,
				vertex_list[(j+2) % vertex_count]->p
			);

			printf("translate([%f,%f,0]) cylinder(r=%f, h=%f, center=true);\n",
				x,
				y,
				hole_rad,
				10.0
			);
		}

		printf("}\n");
	}

	printf("}\n}\n");

	printf("model();\n");


	// For each vertex, extract a small region around the corner
	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];
		const v3_t origin = v->p;

		v3_t avg = {{ 0, 0, 0}};
		find_normal(stl, v, inset_dist, &avg);
		printf("%%translate([%f,%f,%f])", 
			origin.p[0], origin.p[1], origin.p[2]);
		print_normal(&avg);
		printf("rotate([0,180,0]) cylinder(r=15,h=10);\n");

		//avg = v3_norm(avg);

	}

#if 0
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


		// add the screw holes
		make_faces(stl, v, 0, 0, 0, hole_dist, hole_rad, thickness*3);
		printf("} // difference\n");

		printf("}\n");

		refframe_t avg;
		refframe_init(&avg, avg_x, avg_y, avg_z);

		printf("translate([%d,0,12]) render() intersection() {\n", i*30);
		printf("rotate([0,-90,0])");
		print_multmatrix(&avg, 1);
		printf("vertex_%d();\n", i);
		printf("cube([100,100,24], center=true);\n");
		printf("}\n");
		

		//break;
		//if (i == 0) break; // only do one right now
	}

	//printf("translate([0,0,20]) sphere(r=2);\n");
#endif
	return 0;
}
