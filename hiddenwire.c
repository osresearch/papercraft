/** \file
 * Render a hidden wireframe version of an STL file.
 *
 */
// ./hiddenwire --no-hidden --prune 1 -v < nyc-50000.stl  --camera 400,60,-600 --lookat 450,0,-800 --up 0,1,0 --fov 20  > test3.svg

static int debug = 0;

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <err.h>
#include <assert.h>
#include <getopt.h>
#include "v3.h"
#include "tri.h"
#include "camera.h"
#include "svg.h"


static const char usage[] =
"Usage: hiddenwire [options] file.stl > file.svg\n"
"\n"
"Options:\n"
" -h | -? | --help       Help\n"
" -v | --verbose         Enable debugging output\n"
" -c | --camera x,y,z    Camera position\n"
" -l | --lookat x,y,z    Target\n"
" -u | --up x,y,z        Up vector\n"
" -F | --fov deg         Field-of-view angle\n"
" -s | --scale s         Scale factor\n"
" -p | --prune s         Prune lines shorter than s\n"
" --no-backface          Disable backface culling\n"
" --no-coplanar          Disable coplanar merging\n"
" --no-hiddenwire        Disable hidden wire frame removal\n"
"\n";

static const struct option long_options[] = {
	{ "help",			0, NULL, 'h' },
	{ "verbose",			0, NULL, 'v' },
	{ "no-backface",		0, NULL, 'B' },
	{ "no-coplanar",		0, NULL, 'C' },
	{ "no-hiddenwire",		0, NULL, 'H' },
	{ "camera",			1, NULL, 'c' },
	{ "lookat",			1, NULL, 'l' },
	{ "up",				1, NULL, 'u' },
	{ "scale",			1, NULL, 's' },
	{ "prune",			1, NULL, 'p' },
	{ "fov",			1, NULL, 'F' },
	{ NULL,				0, NULL, 0 },
};



#ifndef M_PI
#define 	M_PI   3.1415926535897932384
#endif


typedef struct
{
	char header[80];
	uint32_t num_triangles;
} __attribute__((__packed__))
stl_header_t;


typedef struct
{
	v3_t normal;
	v3_t p[3];
	uint16_t attr;
} __attribute__((__packed__))
stl_face_t;




static inline int
v2_eq(
	const float p0[],
	const float p1[],
	const float eps
)
{
	const float dx = p0[0] - p1[0];
	const float dy = p0[1] - p1[1];

	// are the points within epsilon of each other?
	if (-eps < dx && dx < eps
	&&  -eps < dy && dy < eps)
		return 1;

	// nope, not equal
	return 0;
}





/*
 * Determine if a segment is part of a triangle edge
 */
int
parallel(
	const v3_t * const p00,
	const v3_t * const p01,
	const v3_t * const p10,
	const v3_t * const p11
)
{
	//v3_t v = v3_sub(*p11, *p10);
	v3_t v0 = v3_sub(*p01, *p00);
	v3_t v1 = v3_sub(*p11, *p10);
	v3_t vx = v3_cross(v0, v1);

	float angle = v3_mag2(vx);

	// if the angle is far from zero, definitely not parallel
	if (angle < -EPS || EPS < angle)
		return 0;

	// they might be parallel, figure out if they are the same
	return 1;
}



int v3_parse(v3_t * out, const char * str)
{
	int rc = sscanf(str, "%f,%f,%f",
		&out->p[0],
		&out->p[1],
		&out->p[2]
	);
	if (rc != 3)
		return -1;
	return 0;
}


int onscreen(
	const v3_t * const p,
	const float width,
	const float height
)
{
	if (p->p[0] < -width/2 || width/2 < p->p[0])
		return 0;
	if (p->p[1] < -height/2 || height/2 < p->p[1])
		return 0;
	return 1;
}


int main(
	int argc,
	char ** argv
)
{
	if (argc <= 1)
	{
		fprintf(stderr, "%s", usage);
		return EXIT_FAILURE;
	}

	int opt;
	int do_backface = 1;
	int do_coplanar = 1;
	int do_hidden = 1;
	v3_t eye = { { 100, 0, 0 } };
	v3_t lookat = { { 0, 0, 0 } };
	v3_t up = { { 0, 0, 1 } };
	float scale = 1;
	float fov = 45;
	float prune = 0.1;
	float width = 3000;
	float height = 2000;

	while((opt = getopt_long(argc, argv ,"h?vBCHc:l:s:u:p:F:", long_options, NULL)) != -1)
	{
		switch(opt)
		{
		case 'h' : case '?':
			printf("%s", usage);
			return EXIT_SUCCESS;
		default:
			fprintf(stderr, "%s", usage);
			return EXIT_FAILURE;
		case 'v': debug++; break;

		case 'B': do_backface = 0; break;
		case 'C': do_coplanar = 0; break;
		case 'H': do_hidden = 0; break;

		case 'p': prune = atof(optarg); break;
		case 's': scale = atof(optarg); break;
		case 'F': fov = atof(optarg); break;

		case 'c':
			if (v3_parse(&eye, optarg) < 0)
				return EXIT_FAILURE;
			break;
		case 'l':
			if (v3_parse(&lookat, optarg) < 0)
				return EXIT_FAILURE;
			break;
		case 'u':
			if (v3_parse(&up, optarg) < 0)
				return EXIT_FAILURE;
			break;
		}
	}

	// todo: sanity check fov, scale, etc

	const size_t max_len = 32 << 20;
	uint8_t * const buf = calloc(max_len, 1);
	size_t offset = 0;

	while(1)
	{
		ssize_t rc = read(0, buf+offset, max_len - offset);
		if (rc == -1)
			return EXIT_FAILURE;
		if (rc == 0)
			break;
		offset += rc;
	}

	const stl_header_t * const hdr = (const void*) buf;
	const stl_face_t * const stl_faces = (const void*)(hdr+1);
	const int num_triangles = hdr->num_triangles;

	float coplanar_eps = 0.001;

	if(debug)
	{
		fprintf(stderr, "header: '%s'\n", hdr->header);
		fprintf(stderr, "num: %d\n", num_triangles);
	}


	(void) scale;
	const camera_t * const cam = camera_new(eye, lookat, up, fov);

	printf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0fpx\" height=\"%.0fpx\" viewbox=\"0 0 %.0f %.0f\">\n", width, height, width, height);

	float off_x = width/2;
	float off_y = height/2;
	printf("<g transform=\"translate(%f %f)\">\n", off_x, off_y);

	int rejected = 0;
	tri_t * zlist = NULL;
	seg_t * slist = NULL;
	seg_t * slist_visible = NULL;

	int retained = 0;
	int backface = 0;
	int small_area = 0;
	int behind = 0;
	int offscreen = 0;

	// transform the stl by the camera projection and generate
	// a z-sorted list of triangles
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];

		v3_t s[3];

		for(int j = 0 ; j < 3 ; j++)
		{
			// if any points are behind us, reject
			// this one
			if (!camera_project(cam, &stl->p[j], &s[j]))
			{
				behind++;
				goto reject_early;
			}
		}

		if(debug >= 2)
		for(int j = 0 ; j < 3 ; j++)
		{
			fprintf(stderr, "%+8.1f %+8.1f %+8.1f -> %+8.1f %+8.1f %+8.1f\n",
				stl->p[j].p[0],
				stl->p[j].p[1],
				stl->p[j].p[2],
				s[j].p[0],
				s[j].p[1],
				s[j].p[2]
			);
		}


		tri_t * const tri = tri_new(s, stl->p);

		// reject this face if any of the vertices are behind us
		if (tri->min.p[2] < 0)
		{
			behind++;
			goto reject;
		}

		// do a back-face cull to determine if this triangle
		// is not facing us. we have to determine the orientation
		// from the winding of the new projection
		if (do_backface && tri->normal.p[2] <= 0)
		{
			backface++;
			goto reject;
		}

		// if it has any off-screen coords, reject it
		if (!onscreen(&tri->p[0], width, height)
		||  !onscreen(&tri->p[1], width, height)
		||  !onscreen(&tri->p[2], width, height))
		{
			offscreen++;
			goto reject;
		}

		// prune the small triangles in the screen space
		if (tri_area_2d(tri) < prune)
		{
			small_area++;
			goto reject;
		}

		const float a = v3_dist_2d(&tri->p[0], &tri->p[1]);
		const float b = v3_dist_2d(&tri->p[1], &tri->p[2]);
		const float c = v3_dist_2d(&tri->p[2], &tri->p[0]);
		if( a < prune || b < prune || c < prune)
		{
			small_area++;
			goto reject;
		}

		// it passes the first tests, so insert the triangle
		// into the list and the three line segments
		tri_insert(&zlist, tri);
		retained++;
		continue;

reject:
		tri_delete(tri);
reject_early:
		continue;
	}

	if (debug > 3)
	for(tri_t * t = zlist ; t ; t = t->next)
		tri_print(t);

	if (debug)
		fprintf(stderr, "Retained %d triangles, rejected %d behind, %d offscreen, %d backface, %d small\n", retained, behind, offscreen, backface, small_area);

	// drop any triangles that are totally occluded by another
	// triangle.  this reduces the amount of work for later
	rejected = 0;
	for(tri_t * t = zlist ; t ; t = t->next)
	{
		tri_t * t2_next;
		for(tri_t * t2 = zlist ; t2 ; t2 = t2_next)
		{
			t2_next = t2->next;
			if (t == t2)
				continue;

			if (!tri_behind(t, t2))
				continue;

			// t2 is occluded by t, remove it from the list
			rejected++;
			tri_delete(t2);
		}
		
	}
	if (debug)
		fprintf(stderr, "Rejected %d fully occluded triangles\n", rejected);


	// generate a list of segments, dropping any coplanar ones
	rejected = 0;
	for(tri_t * t = zlist ; t ; t = t->next)
	{
		unsigned matches = 0;

		if(do_coplanar)
		for(tri_t * t2 = zlist ; t2 ; t2 = t2->next)
		{
			if (t == t2)
				continue;

			const int edge = tri_coplanar(t, t2, coplanar_eps);
			if (edge < 0)
				continue;
			matches |= 1 << edge;
		}
		
		for(int j = 0 ; j < 3 ; j++)
		{
			// drop any that are coplanar
			if (matches & (1 << j))
			{
				rejected++;
				continue;
			}

			seg_t * s = seg_new(t->p[j], t->p[(j+1) % 3]);
			s->next = slist;
			slist = s;
		}
	}

	if (debug)
		fprintf(stderr, "Rejected %d coplanar segments\n", rejected);


	// we now have a z-sorted list of triangles
	rejected = 0;

	// compute how many we actuall have remaining
	int remaining = 0;
	if (debug)
	{
		for(seg_t * s = slist ; s ; s = s->next)
			remaining++;
		fprintf(stderr, "%d segments remain to process\n", remaining);
	}


	if(do_hidden)
	{
		// work on each segment, intersecting it with all of the triangles
		int processed = 0;
		while(slist)
		{
			if (debug && ++processed % 1000 == 0)
				fprintf(stderr, "Hidden %d\n", processed);

			seg_t * s = slist;
			slist = s->next;

			tri_seg_hidden(zlist, s, &slist_visible);
		}
	} else {
		// don't do any intersection tests
		slist_visible = slist;
		slist = NULL;
	}

	// display all of the visible segments
	for(seg_t * s = slist_visible ; s ; s = s->next)
	{
		svg_line("#FF0000", s->p[0].p, s->p[1].p, 1);
	}


	printf("</g>\n");
	printf("</svg>\n");

	return 0;
}
