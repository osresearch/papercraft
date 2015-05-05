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
#include <getopt.h>
#include "v3.h"
#include "stl_3d.h"

static FILE * output;
static int verbose;


static void
print_multmatrix(
	const refframe_t * const ref,
	const int transpose
)
{
	fprintf(output, "multmatrix(m=["
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
	const v3_t * normal,
	int show_model
)
{
	const float x = normal->p[0];
	const float y = normal->p[1];
	const float z = normal->p[2];
	const double length = sqrt(x*x+y*y+z*z);
	const double b = acos(z / length);
	const double c = x == 0 ? sign(y)*90 : atan2(y,x);

	if (!show_model)
	{
		fprintf(output, "rotate([0,%f,0])", -b*180/M_PI);
		fprintf(output, "rotate([0,0,%f])", -c*180/M_PI);
	} else {
		fprintf(output, "rotate([%f,%f,%f])\n", 0.0, b * 180 / M_PI, c * 180 / M_PI);
	}
}


static void
find_normal(
	const stl_3d_t * const stl,
	const stl_vertex_t * const v,
	const float inset_distance,
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
			refframe_inset(&ref, inset_distance, &x, &y, p1, p2, p3);

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

			//*avg = v3_add(*avg, ref.z);
			*avg = v3_add(*avg, v3_norm(v3_sub(ref.origin, hole)));
			//*avg = v3_add(*avg, (v3_sub(hole, ref.origin)));
		}
	}

	free(face_used);
	free(vertex_list);
}


static struct option long_options[] = 
{
	{ "verbose",	no_argument,	   0, 'v' },
	{ "model",	no_argument,	   0, 'm' },
	{ "inset",	required_argument, 0, 'i' },
	{ "radius",     required_argument, 0, 'r' },
	{ "input",      required_argument, 0, 'I' },
	{ "output",	required_argument, 0, 'O' },
	{ 0, 0, 0, 0 },
};


static void
usage(
	FILE * const out
)
{
	fprintf(out,
		"Usage: corners [options] -I stl-binary.stl > corners.scad\n"
		"Options:\n"
		" -v | --verbose      Enable verbosity\n"
		" -i | --inset N      Inset mm\n"
		" -r | --radius N     Hole radius mm\n"
		" -I | --input file   Read binary STL from file\n"
		" -O | --output file  Write SVG to file\n"
		" -m | --model        Generate a 3D model instead of corners\n"
		"\n"
	);
}



int
main(
	int argc,
	char ** argv
)
{
	double inset_distance = 5;
	double hole_radius = 1.15;
	const char * input_file = NULL;
	const char * output_file = NULL;
	int show_model = 0;
	int option_index = 0;

	while (1)
	{
		const int c = getopt_long(
			argc,
			argv,
			"vmI:r:i:O:",
			long_options,
			&option_index
		);
		if (c == -1)
			break;
		switch(c)
		{
		case 'm': show_model = 1; break;
		case 'v': verbose++; break;
		case 'i': inset_distance = atof(optarg); break;
		case 'r': hole_radius = atof(optarg); break;
		case 'I': input_file = optarg; break;
		case 'O': output_file = optarg; break;
		case 'h': case '?':
			usage(stdout);
			return 0;
		default:
			usage(stderr);
			return -1;
		}
	}

	int input_fd;
	if (!input_file)
	{
		fprintf(stderr, "Input STL must be specified\n");
		return -1;
	} else {
		input_fd = open(input_file, O_RDONLY);
		if (input_fd < 0)
		{
			perror(input_file);
			return -1;
		}
	}

	if (!output_file)
	{
		output_file = "stdout";
		output = stdout;
	} else {
		output = fopen(output_file, "w");
		if (!output)
		{
			perror(output_file);
			return -1;
		}
	}

	stl_3d_t * const stl = stl_3d_parse(input_fd);
	if (!stl)
	{
		fprintf(stderr, "%s: Unable to parse STL\n", input_file);
		return EXIT_FAILURE;
	}
	close(input_fd);


	if (verbose)
		fprintf(stderr,
			"%s: %d faces, %d vertex\n",
			input_file,
			stl->num_face,
			stl->num_vertex
		);

	fprintf(output, "module model() {\n"
		"render() difference() {\n"
	 	"import(\"%s\");\n",
		input_file
	);
	//printf("%%model();\n");

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

		fprintf(output, "translate([%f,%f,%f])",
			f->vertex[0]->p.p[0],
			f->vertex[0]->p.p[1],
			f->vertex[0]->p.p[2]
		);
			
		print_multmatrix(&ref, 0);
		fprintf(output, "{\n");

		// generate a bolt hole for each non-copolanar corner
		for (int j = 0 ; j < vertex_count ; j++)
		{
			double x, y;
			refframe_inset(
				&ref,
				inset_distance,
				&x,
				&y,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p,
				vertex_list[(j+2) % vertex_count]->p
			);

			fprintf(output, "translate([%f,%f,0]) cylinder(r=%f, h=%f, center=true);\n",
				x,
				y,
				hole_radius,
				10.0
			);
		}

		fprintf(output, "}\n");
	}

	fprintf(output, "}\n}\n");

	if (show_model)
		fprintf(output, "model();\n");


	// For each vertex, extract a small region around the corner
	const int div = sqrt(stl->num_vertex);
	const double spacing = 32;

	for(int i = 0 ; i < stl->num_vertex ; i++)
	{
		const stl_vertex_t * const v = &stl->vertex[i];
		const v3_t origin = v->p;

		v3_t avg = {{ 0, 0, 0}};
		find_normal(stl, v, inset_distance, &avg);

		if (!show_model)
		{
			fprintf(output, "translate([%f,%f,20])", (i/div)*spacing, (i%div)*spacing);
			fprintf(output, "render() intersection()");
		}

		fprintf(output, "{\n");

		//printf("%%\n");
		if (!show_model)
		{
			print_normal(&avg, show_model);
			fprintf(output, "translate([%f,%f,%f])", 
				-origin.p[0], -origin.p[1], -origin.p[2]);
			fprintf(output, "model();\n");
			fprintf(output, "translate([0,0,-20]) cylinder(r=15,h=20);\n");
		} else {
			fprintf(output, "translate([%f,%f,%f])", 
				origin.p[0], origin.p[1], origin.p[2]);
			print_normal(&avg, show_model);
			fprintf(output, "%%translate([0,0,-20]) cylinder(r=15,h=20);\n");
		}

		//avg = v3_norm(avg);

		fprintf(output, "}\n");
	}

	return 0;
}
