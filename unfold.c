/** \file
 * Unfold an STL file into a set of laser-cutable polygons.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>

#ifndef M_PI
#define 	M_PI   3.1415926535897932384
#endif

#define EPS 0.0001

typedef struct
{
	char header[80];
	uint32_t num_triangles;
} __attribute__((__packed__))
stl_header_t;

typedef struct
{
	float p[3];
} v3_t;

typedef struct
{
	v3_t normal;
	v3_t p[3];
	uint16_t attr;
} __attribute__((__packed__))
stl_face_t;


typedef struct face face_t;

struct face
{
	float sides[3];
	face_t * next[3];
	int next_edge[3];
	int coplanar[3];
	int used;
};


static int
v3_eq(
	const v3_t * v1,
	const v3_t * v2
)
{
	float dx = v1->p[0] - v2->p[0];
	float dy = v1->p[1] - v2->p[1];
	float dz = v1->p[2] - v2->p[2];

	if (-EPS < dx && dx < EPS
	&&  -EPS < dy && dy < EPS
	&&  -EPS < dz && dz < EPS)
		return 1;

	return 0;
}


static int
edge_eq(
	const stl_face_t * const t0,
	const stl_face_t * const t1,
	int e0,
	int e1
)
{
	const v3_t * const v0 = &t0->p[e0];
	const v3_t * const v1 = &t0->p[e1];

	if (v3_eq(v0, &t1->p[0]) && v3_eq(v1, &t1->p[1]))
		return 1;
	if (v3_eq(v0, &t1->p[1]) && v3_eq(v1, &t1->p[0]))
		return 1;

	if (v3_eq(v0, &t1->p[0]) && v3_eq(v1, &t1->p[2]))
		return 1;
	if (v3_eq(v0, &t1->p[2]) && v3_eq(v1, &t1->p[0]))
		return 1;

	if (v3_eq(v0, &t1->p[1]) && v3_eq(v1, &t1->p[2]))
		return 1;
	if (v3_eq(v0, &t1->p[2]) && v3_eq(v1, &t1->p[1]))
		return 1;

	return 0;
}


/* Compare two edges in two triangles.
 *
 * note that if the windings are all the same, the edges will
 * compare in the opposite order (for example, the edge from 0 to 1
 * compares to the edge from 2 to 1 in the other triangle).
 */
static int
edge_eq2(
	const stl_face_t * const t0,
	const stl_face_t * const t1,
	int e0,
	int e1
)
{
	const v3_t * const v00 = &t0->p[e0];
	const v3_t * const v01 = &t0->p[(e0+1) % 3];

	const v3_t * const v10 = &t1->p[e1];
	const v3_t * const v11 = &t1->p[(e1+1) % 3];

	if (v3_eq(v00, v11) && v3_eq(v01, v10))
		return 1;

	return 0;
}


double
v3_len(
	const v3_t * const v0,
	const v3_t * const v1
)
{
	float dx = v0->p[0] - v1->p[0];
	float dy = v0->p[1] - v1->p[1];
	float dz = v0->p[2] - v1->p[2];

	return sqrt(dx*dx + dy*dy + dz*dz);
}


/** recursively try to fix up the triangles.
 *
 * returns 0 if this should be unwound, 1 if was successful
 */
int
recurse(
	face_t * const f,
	int start_edge
)
{
	static int depth;

	depth++;

	// flag that we are looking into this one
	f->used = 1;

	// print out a svg group for this triangle, starting with
	// the incoming edge
	float a = f->sides[(start_edge + 0) % 3];
	float c = f->sides[(start_edge + 1) % 3];
	float b = f->sides[(start_edge + 2) % 3];
	float x2 = (a*a + b*b - c*c) / (2*a);
	float y2 = sqrt(b*b - x2*x2);

	printf("<polyline points=\"0,0 %f,0 %f,%f\ 0,0\" fill=\"none\" stroke=\"#FF0000\" />\n",
		a,
		x2,
		y2
	);
	//printf("%p %d %f %f %f\n", f, start_edge, f->sides[0], f->sides[1], f->sides[2]);

	// for each edge, find the triangle that matches
	for (int edge = 0 ; edge < 3 ; edge++)
	{
		face_t * const f2 = f->next[(edge+start_edge) % 3];
		if (f2->used)
			continue;

		// create a group that translates and rotates
		// such that it lines up with this edge
		float trans_x, trans_y, rotate;
		if (edge == 0)
		{
			trans_x = a;
			trans_y = 0;
			rotate = 180;
		} else
		if (edge == 1)
		{
			trans_x = x2;
			trans_y = y2;
			rotate = -atan2(y2, a-x2) * 180 / M_PI;
		} else
		if (edge == 2)
		{
			trans_x = 0;
			trans_y = 0;
			rotate = atan2(y2, x2) * 180 / M_PI;
		}

		printf("<!-- edge %d --><g transform=\"translate(%f,%f) rotate(%f)\">\n",
			edge,
			trans_x,
			trans_y,
			rotate
		);

		recurse(f2, f->next_edge[(edge+start_edge) % 3]);

		printf("</g>\n");
	}

	// no success
	return 0;
}


int
coplanar_check(
	const stl_face_t * const f1,
	const stl_face_t * const f2
)
{
	// no, for now
	return 0;
}



int main(void)
{
	const size_t max_len = 1 << 20;
	uint8_t * const buf = calloc(max_len, 1);

	ssize_t rc = read(0, buf, max_len);
	if (rc == -1)
		return EXIT_FAILURE;

	const stl_header_t * const hdr = (const void*) buf;
	const stl_face_t * const stl_faces = (const void*)(hdr+1);
	const int num_triangles = hdr->num_triangles;

	fprintf(stderr, "header: '%s'\n", hdr->header);
	fprintf(stderr, "num: %d\n", num_triangles);

	face_t * const faces = calloc(num_triangles, sizeof(*faces));

	// convert the stl triangles into faces
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];
		face_t * const f = &faces[i];

		f->sides[0] = v3_len(&stl->p[0], &stl->p[1]);
		f->sides[1] = v3_len(&stl->p[1], &stl->p[2]);
		f->sides[2] = v3_len(&stl->p[2], &stl->p[0]);
	}

	// look to see if there is a matching point
	// in the faces that we've already built
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];
		face_t * const f = &faces[i];

		for (int j = 0 ; j < num_triangles ; j++)
		{
			if (i == j)
				continue;

			const stl_face_t * const stl2 = &stl_faces[j];
			face_t * const f2 = &faces[j];

			for (int edge = 0 ; edge < 3 ; edge++)
			{
				if (f->next[edge])
					continue;

				for (int edge2 = 0 ; edge2 < 3 ; edge2++)
				{
					if (f2->next[edge2])
						continue;

					if (!edge_eq2(stl, stl2, edge, edge2))
						continue;

					f->next[edge] = f2;
					f->next_edge[edge] = edge2;
					f2->next[edge2] = f;
					f2->next_edge[edge2] = edge;

					f->coplanar[edge] =
					f2->coplanar[edge2] = coplanar_check(stl, stl2);
				}
			}
		}

		// all three edges should be matched
		if (f->next[0] && f->next[1] && f->next[2])
			continue;
		fprintf(stderr, "%d: missing edges?\n", i);
	}

	// we now have a graph that shows the connection between
	// all of the faces and their sizes. start converting them

	printf("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");
	//for (int i = 0 ; i < num_triangles ; i++)
		recurse(&faces[0], 0);

	printf("</svg>\n");

	return 0;
}
