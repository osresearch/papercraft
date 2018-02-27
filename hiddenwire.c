/** \file
 * Render a hidden wireframe version of an STL file.
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
#include "v3.h"
#include "camera.h"


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

static int debug = 0;

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


typedef struct _tri_t tri_t;
struct _tri_t
{
	v3_t p[3]; // camera space
	v3_t normal; // camera space
	v3_t normal_xyz; // original xyz space
	float min[3]; // camera space
	float max[3]; // camera space
	tri_t * next;
	tri_t ** prev;
};


typedef struct _seg_t seg_t;
struct _seg_t {
	v3_t p[2];
	v3_t src[2];
	seg_t * next;
};



void
svg_line(
	const char * color,
	const float * p1,
	const float * p2,
	float thick
)
{
	printf("<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%s\" stroke-width=\"%.1fpx\"/>\n",
		p1[0],
		p1[1],
		p2[0],
		p2[1],
		color,
		thick
	);
}


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


/** Compute the points of intersection for two segments in 2d, and z points.
 *
 * This is a specialized ray intersection algorithm for the
 * hidden wire-frame removal code that computes the intersection
 * points for two rays (in 2D, "orthographic") and then computes
 * the Z depth for the intersections along each of the segments.
 *
 * Returns -1 for non-intersecting, otherwise a ratio of how far
 * along the intersection is on the l0.
 */
float
hidden_intersect(
	const v3_t * const p0,
	const v3_t * const p1,
	const v3_t * const p2,
	const v3_t * const p3,
	v3_t * const l0_int,
	v3_t * const l1_int
)
{
	const float p0_x = p0->p[0];
	const float p0_y = p0->p[1];
	const float p0_z = p0->p[2];
	const float p1_x = p1->p[0];
	const float p1_y = p1->p[1]; 
	const float p1_z = p1->p[2];
	const float p2_x = p2->p[0];
	const float p2_y = p2->p[1];
	const float p2_z = p2->p[2];
	const float p3_x = p3->p[0];
	const float p3_y = p3->p[1];
	const float p3_z = p3->p[2];

	const float s1_x = p1_x - p0_x;
	const float s1_y = p1_y - p0_y;
	const float s2_x = p3_x - p2_x;
	const float s2_y = p3_y - p2_y;

	// compute r x s
	const float d = -s2_x * s1_y + s1_x * s2_y;

	// if they are close to parallel, then we do not need to check
	// for intersection (we define that as "non-intersecting")
	if (-EPS < d && d < EPS)
		return -1;

	// Compute how far along each line they would interesect
	const float r0 = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / d;
	const float r1 = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / d;

	// if they are not within the ratio (0,1), then the intersecton occurs
	// outside of the segments and is not of concern
	if (r0 < 0 || r0 > 1)
		return -1;
	if (r1 < 0 || r1 > 1)
		return -1;

	// Collision detected with the segments
if(0) fprintf(stderr, "collision: %.0f,%.0f,%.0f->%.0f,%.0f,%.0f %.0f,%.0f,%.0f->%.0f,%.0f,%.0f == %.3f,%.3f\n",
		p0_x, p0_y, p0_z,
		p1_x, p1_y, p1_z,
		p2_x, p2_y, p2_z,
		p3_x, p3_y, p2_z,
		r0,
		r1
	);

	const float ix = p0_x + (r0 * s1_x);
	const float iy = p0_y + (r0 * s1_y);

	// compute the z intercept for each on the two different coordinates
	if(l0_int)
	{
		*l0_int = (v3_t){{
			ix,
			iy,
			p0_z + r0 * (p1_z - p0_z)
		}};
	}

	if(l1_int)
	{
		*l1_int = (v3_t){{
			ix,
			iy,
			p2_z + r1 * (p3_z - p2_z)
		}};
	}

	return r0;
}





tri_t *
tri_new(
	const v3_t * p_cam,
	const v3_t * p_xyz
)
{
	tri_t * const t = calloc(1, sizeof(*t));
	if (!t)
		return NULL;
	for(int i = 0 ; i < 3  ; i++)
		t->p[i] = p_cam[i];

	// precompute the normals
	t->normal = v3_norm(v3_cross(
		v3_sub(t->p[1], t->p[0]),
		v3_sub(t->p[2], t->p[1])
	));
	t->normal_xyz = v3_norm(v3_cross(
		v3_sub(p_xyz[1], p_xyz[0]),
		v3_sub(p_xyz[2], p_xyz[1])
	));


	// compute the bounding box for the triangle in camera space
	for(int j = 0 ; j < 3 ; j++)
	{
		t->min[j] = min(min(t->p[0].p[j], t->p[1].p[j]), t->p[2].p[j]);
		t->max[j] = max(max(t->p[0].p[j], t->p[1].p[j]), t->p[2].p[j]);
	}

	return t;
}


// insert a triangle into our z-sorted list
void
tri_insert(
	tri_t ** zlist,
	tri_t * t
)
{
	while(1)
	{
		tri_t * const iter = *zlist;
		if (!iter)
			break;

		// check to see if our new triangle is closer than
		// the current triangle
		if(iter->min[2] > t->min[2])
			break;

		zlist = &(iter->next);
	}

	// either we reached the end of the list,
	// or we have found where our new triangle is sorted
	t->next = *zlist;
	*zlist = t;
	if (t->next)
		t->next->prev = &t->next;
}


void
tri_delete(tri_t * t)
{
	if (t->next)
		t->next->prev = t->prev;
	if (t->prev)
		*(t->prev) = t->next;

	t->next = NULL;
	t->prev = NULL;
	free(t);
}


seg_t *
seg_new(
	const v3_t p0,
	const v3_t p1
)
{
	seg_t * const s = calloc(1, sizeof(*s));
	if (!s)
		return NULL;
	s->p[0] = p0;
	s->p[1] = p1;
	s->src[0] = p0;
	s->src[1] = p1;
	s->next = NULL;

	return s;
}


void
seg_print(
	const seg_t * const s
)
{
	fprintf(stderr, "%.0f,%.0f -> %.0f,%.0f (was %.0f,%.0f -> %.0f,%.0f\n",
		s->p[0].p[0],
		s->p[0].p[1],
		s->p[1].p[0],
		s->p[1].p[1],
		s->src[0].p[0],
		s->src[0].p[1],
		s->src[1].p[0],
		s->src[1].p[1]
	);
}

void
tri_print(
	const tri_t * const t
)
{
	fprintf(stderr, "%.0f,%.0f,%.0f %.0f,%.0f,%.0f %.0f,%.0f,%.0f norm %.3f,%.3f,%.3f\n",
		t->p[0].p[0],
		t->p[0].p[1],
		t->p[0].p[2],
		t->p[1].p[0],
		t->p[1].p[1],
		t->p[1].p[2],
		t->p[2].p[0],
		t->p[2].p[1],
		t->p[2].p[2],
		t->normal.p[0],
		t->normal.p[1],
		t->normal.p[2]
	);
}


/* Check if two triangles are coplanar and share an edge.
 *
 * Returns -1 if not coplanar, 0-2 for the edge in t0 that they share.
 */
int
tri_coplanar(
	const tri_t * const t0,
	const tri_t * const t1,
	const float coplanar_eps
)
{
	// the two normals must be parallel-enough
	const float angle = v3_mag(v3_sub(t0->normal_xyz, t1->normal_xyz));
	if (angle < -coplanar_eps || +coplanar_eps < angle)
		return -1;
	
	// find if there are two points shared
	unsigned matches = 0;
	for(int i = 0 ; i < 3 ; i++)
	{
		for(int j = 0 ; j < 3 ; j++)
		{
			if (!v3_eq(&t0->p[i], &t1->p[j]))
				continue;
			matches |= 1 << i;
			break;
		}
	}

	switch(matches)
	{
	case 0x3: return 0;
	case 0x6: return 1;
	case 0x5: return 2;
	case 0x7:
		fprintf(stderr, "uh, three points match?\n");
		tri_print(t0);
		tri_print(t1);
		return -1;
	default:
		// no shared edge
		return -1;
	}
}


/*
 * Find the Z point of an XY coordinate in a triangle.
 *
 * p can be written as a combination of t01 and t02,
 * p - t0 = a * (t1 - t0) + b * (t2 - t0)
 * setting t0 to 0, this becomes:
 * p = a * t1 + b * t2
 * which is two equations with two unknowns
 */
int
tri_find_z(
	const tri_t * const t,
	const v3_t * const p,
	float * const zout
)
{
	const float t1x = t->p[1].p[0] - t->p[0].p[0];
	const float t1y = t->p[1].p[1] - t->p[0].p[1];
	const float t1z = t->p[1].p[2] - t->p[0].p[2];
	const float t2x = t->p[2].p[0] - t->p[0].p[0];
	const float t2y = t->p[2].p[1] - t->p[0].p[1];
	const float t2z = t->p[2].p[2] - t->p[0].p[2];
	const float px = p->p[0] - t->p[0].p[0];
	const float py = p->p[1] - t->p[0].p[1];

	const float a = (px * t2y - py * t2x) / (t1x * t2y - t2x * t1y);
	const float b = (px * t1y - py * t1x) / (t2x * t1y - t1x * t2y);

	const float z = t->p[0].p[2] + a * t1z + b * t2z;

	if (zout)
		*zout = z;

	return 0 <= a && 0 <= b && a + b <= 1;
}


/*
 * Recursive algorithm:
 * Given a line segment and a list of triangles,
 * find if the line segment crosses any triangle.
 * If it crosses a triangle the segment will be shortened
 * and an additional one might be created.
 * Recusively try intersecting the new segment (starting at the same triangle)
 * and then continue trying the shortened segment.
 */

void
tri_seg_intersect(
	const tri_t * zlist,
	seg_t * s,
	seg_t ** slist_visible
)
{
	const float p0z = s->p[0].p[2];
	const float p1z = s->p[1].p[2];
	const float seg_max_z = max(p0z, p1z);

	// avoid processing empty segments
	const float seg_len = v3_len(&s->p[0], &s->p[1]);
	if (seg_len < EPS)
		return;

static int recursive;
recursive++;
//fprintf(stderr, "%d: processing segment ", recursive++); seg_print(s);

	for( const tri_t * t = zlist ; t ; t = t->next )
	{
		// if the segment is closer than the triangle,
		// then we no longer have to check any further into
		// the zlist (it is sorted by depth).
		if (seg_max_z <= t->min[2])
			break;

		// make sure that we're not comparing to our own triangle
		// or one that shares an edge with us (which might be in
		// a different order)
		if (v2_eq(s->src[0].p, t->p[0].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[1].p, 0.0005))
			continue;
		if (v2_eq(s->src[0].p, t->p[1].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[2].p, 0.0005))
			continue;
		if (v2_eq(s->src[0].p, t->p[2].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[0].p, 0.0005))
			continue;
		if (v2_eq(s->src[0].p, t->p[1].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[0].p, 0.0005))
			continue;
		if (v2_eq(s->src[0].p, t->p[2].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[1].p, 0.0005))
			continue;
		if (v2_eq(s->src[0].p, t->p[0].p, 0.0005)
		&&  v2_eq(s->src[1].p, t->p[2].p, 0.0005))
			continue;

		float z0, z1;
		int inside0 = tri_find_z(t, &s->p[0], &z0);
		int inside1 = tri_find_z(t, &s->p[1], &z1);

		// if both are inside but the segment is infront of the
		// triangle, then we retain the segment.
		// otherwies we discard the segment
		if (inside0 && inside1)
		{
			if (s->p[0].p[2] <= z0
			&&  s->p[1].p[2] <= z1)
				continue;
			recursive--;
			return;
		}

		// split the segment for each intersection with the
		// triangle segments and add it to the work queue.
		int intersections = 0;
		v3_t is[3] = {}; // 3d point of segment intercept
		v3_t it[3] = {}; // 3d point of triangle intercept

		for(int j = 0 ; j < 3 ; j++)
		{
			float ratio = hidden_intersect(
				&s->p[0], &s->p[1],
				&t->p[j], &t->p[(j+1)%3],
				&is[intersections],
				&it[intersections]
			);

			if (ratio < 0)
				continue;

			intersections++;
		}

		// if none of them intersect, we keep looking
		if (intersections == 0)
			continue;

		if (intersections == 3)
		{
			// this likely means that the triangle is very, very
			// small, so let's just throw away this line segment
			recursive--;
			return;
		}


		if (intersections == 2)
		{
			// figure out how far it is to each of the intersections
			const float d00 = v3_len(&s->p[0], &is[0]);
			const float d01 = v3_len(&s->p[0], &is[1]);
			const float d10 = v3_len(&s->p[1], &is[0]);
			const float d11 = v3_len(&s->p[1], &is[1]);

			// discard segments that have two interesections that match
			// the segment exactly (distance from segment ends to
			// intersection point close enough to zero).
			if (d00 < EPS && d11 < EPS)
			{
				recursive--;
				return;
			}
			if (d01 < EPS && d10 < EPS)
			{
				recursive--;
				return;
			}

			// if the segment intersection is closer than the triangle,
			// then we do nothing. degenerate cases are not handled
			if (d00 <= d01
			&& is[0].p[2] <= it[0].p[2]
			&& is[1].p[2] <= it[1].p[2])
				continue;
			if (d00 > d01
			&& is[1].p[2] <= it[0].p[2]
			&& is[0].p[2] <= it[1].p[2])
				continue;

			// segment is behind the triangle,
			// we have to create a new segment
			// and shorten the existing segment
			// find the two intersections that we have
			// update the src field

			// we need to create a new segment
			seg_t * news;
			if (d00 < d01)
			{
				// split from p0 to ix0
				news = seg_new(s->p[0], is[0]);
				news->src[0] = s->src[0];
				news->src[1] = s->src[1];
				s->p[0] = is[1];
			} else {
				// split from p0 to ix1
				news = seg_new(s->p[0], is[1]);
				news->src[0] = s->src[0];
				news->src[1] = s->src[1];
				s->p[0] = is[0];
			}

			// recursively start splitting the new segment
			// starting at the next triangle down the z-depth
			tri_seg_intersect(zlist->next, news, slist_visible);

			// continue splitting our current segment
			continue;
		}

		if (intersections == 1)
		{
			// if there is an intersection, but the segment intercept
			// is closer than the triangle intercept, then no problem.
			// we do not bother with degenerate cases of intersecting
			// triangles
			if (is[0].p[2] <= it[0].p[2]
			&&  is[1].p[2] <= it[0].p[2])
			{
				//svg_line("#00FF00", s->p[0].p, s->p[1].p, 10);
				continue;
			}

			if (inside0)
			{
				// shorten it on the 0 side
				s->p[0] = is[0];
// huh? shouldn't we process this one?
return;
				continue;
			} else
			if (inside1)
			{
				// shorten it on the 1 side
				s->p[1] = is[0];
// huh? shouldn't we process this one?
return;
				continue;
			} else {
				// both outside, but an intersection?
				// split at that point and hope for the best
				seg_t * const news = seg_new(s->p[0], is[0]);
				news->src[0] = s->src[0];
				news->src[1] = s->src[1];
				s->p[0] = is[0];

				tri_seg_intersect(zlist->next, news, slist_visible);
				// continue splitting our current segment
				continue;
			}
		}
	}

	// if we've reached here the segment is visible
	// and should be added to the visible list
	s->next = *slist_visible;
	*slist_visible = s;
	recursive--;
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


	const camera_t * const cam = camera_new(eye, lookat, up, fov, scale);

	printf("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");

	float off_x = 0;
	float off_y = 0;
	printf("<g transform=\"translate(%f %f)\">\n", off_x, off_y);

	int rejected = 0;
	tri_t * zlist = NULL;
	seg_t * slist = NULL;
	seg_t * slist_visible = NULL;

	int retained = 0;

	// transform the stl by the camera projection and generate
	// a z-sorted list of triangles
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];

		v3_t s[3];

		for(int j = 0 ; j < 3 ; j++)
			camera_project(cam, &stl->p[j], &s[j]);

		if(debug >= 2)
		fprintf(stderr, "%.3f,%.3f,%.3f -> %.1f,%.1f,%.1f\n",
			stl->p[0].p[0],
			stl->p[0].p[1],
			stl->p[0].p[2],
			s[0].p[0],
			s[0].p[1],
			s[0].p[2]
		);

		tri_t * const tri = tri_new(s, stl->p);

		// reject this face if any of the vertices are behind us
		if (tri->min[2] < 0)
			goto reject;

		// do a back-face cull to determine if this triangle
		// is not facing us. we have to determine the orientation
		// from the winding of the new projection
		if (do_backface && tri->normal.p[2] <= 0)
			goto reject;

		retained++;

		// it passes the first tests, so insert the triangle
		// into the list and the three line segments
		tri_insert(&zlist, tri);
		continue;

reject:
		tri_delete(tri);
		rejected++;
	}

	if (debug)
		fprintf(stderr, "Retained %d, rejected %d triangles\n", retained, rejected);


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

	if(do_hidden)
	{
		// work on each segment, intersecting it with all of the triangles
		int processed = 0;
		while(slist)
		{
			if (++processed % 100 == 0)
				fprintf(stderr, "Hidden %d\n", processed);

			seg_t * s = slist;
			slist = s->next;

			tri_seg_intersect(zlist, s, &slist_visible);
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

	if (debug)
		fprintf(stderr, "Occluded %d triangles\n", rejected);


	printf("</g>\n");
	printf("</svg>\n");

	return 0;
}
