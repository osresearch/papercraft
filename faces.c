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
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "v3.h"
#include "stl_3d.h"


static FILE * output;
static int verbose;

static const char * stroke_string
	= "stroke-width=\"1px\" fill=\"none\"";

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

	fprintf(output, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%s\" %s/>\n",
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
	fprintf(output, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"%s\" %s/>\n",
		x,
		y,
		rad,
		color,
		stroke_string
	);
}


static struct option long_options[] = 
{
	{ "verbose",	no_argument,	   0, 'v' },
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
		"Usage: faces [options] < stl-binary.stl > faces.svg\n"
		"Options:\n"
		" -v | --verbose      Enable verbosity\n"
		" -i | --inset N      Inset mm\n"
		" -r | --radius N     Hole radius mm\n"
		" -I | --input file   Read binary STL from file\n"
		" -O | --output file  Write SVG to file\n"
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
	int option_index = 0;

	while (1)
	{
		const int c = getopt_long(
			argc,
			argv,
			"vI:r:i:O:",
			long_options,
			&option_index
		);
		if (c == -1)
			break;
		switch(c)
		{
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
	 	input_fd= STDIN_FILENO;
		input_file = "stdin";
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

	int * const face_used
		= calloc(sizeof(*face_used), stl->num_face);

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	fprintf(output, "<svg xmlns=\"http://www.w3.org/2000/svg\">\n");
	fprintf(output, "<g transform=\"scale(3.543307)\"><!-- scale to mm -->\n");

	const stl_vertex_t ** const vertex_list
		= calloc(sizeof(*vertex_list), stl->num_vertex);

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

		fprintf(stderr, "%s: %d: %d vertices\n", input_file, i, vertex_count);

		// generate a refernce frame based on this face
		refframe_t ref;
		refframe_init(&ref,
			f->vertex[0]->p,
			f->vertex[1]->p,
			f->vertex[2]->p
		);
		fprintf(output, "<!-- face %d --><g>\n", i);

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

		fprintf(output, "</g>\n");
	}

	fprintf(output, "</g></svg>\n");

	fclose(output);

	return 0;
}
