/** \file
 * Generate an OpenSCAD with cubes for each edge
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

#ifndef M_PI
#define 	M_PI   3.1415926535897932384
#endif

static int debug = 0;
static int draw_labels = 0;

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


#define MAX_VERTEX 32

typedef struct stl_vertex stl_vertex_t;

struct stl_vertex
{
	v3_t p;
	int num_edges;
	stl_vertex_t * edges[MAX_VERTEX];
};



/* Returns 1 for ever edge in f1 that is shared with f2.
 */
int
coplanar_check(
	const stl_face_t * const f1,
	const stl_face_t * const f2
)
{
	// Verify that there are three matching points
	int match[3] = {0,0,0};
	for (int i = 0 ; i < 3 ; i++)
	{
		for (int j = 0 ; j < 3 ; j++)
			if (v3_eq(&f1->p[i], &f2->p[j]))
				match[i] = 1;
	}

	uint8_t mask = 0;

	if (match[0] && match[1])
		mask = 1;
	if (match[1] && match[2])
		mask = 2;
	if (match[2] && match[0])
		mask = 4;

	// otherwise they do not share enough points
	if (mask == 0)
		return 0;

	// find the four distinct points
	v3_t x1 = f1->p[0];
	v3_t x2 = f1->p[1];
	v3_t x3 = f1->p[2];
	v3_t x4;

	for (int i = 0 ; i < 3 ; i++)
	{
		x4 = f2->p[i];
		if (v3_eq(&x1, &x4))
			continue;
		if (v3_eq(&x2, &x4))
			continue;
		if (v3_eq(&x3, &x4))
			continue;
		break;
	}

	// (x3-x1) . ((x2-x1) X (x4-x3)) == 0
	v3_t dx31 = v3_sub(x3, x1);
	v3_t dx21 = v3_sub(x2, x1);
	v3_t dx43 = v3_sub(x4, x3);
	v3_t cross = v3_cross(dx21, dx43);
	float dot = v3_dot(dx31, cross);
	
	int check = -EPS < dot && dot < +EPS;

	// if the dot product is not close enough to zero, they
	// are not coplanar.
	if (!check)
		return 0;

	// coplanar! return the shared edge mask
	return mask;
}


static inline float
sign(
	const float x
)
{
	if (x < 0)
		return -1;
	if (x > 0)
		return +1;
	return 0;
}



/**
 * Add a vector to the list of edges if it is not already present
 * and if it is not coplanar with other ones.
 * Note that if it is coplanar, but "outside" the other edges then it
 * will replace the inside one.
 */
void
stl_edge_insert(
	stl_vertex_t * const v1,
	stl_vertex_t * const v2
)
{
	for (int i = 0 ; i < v1->num_edges ; i++)
	{
		const v3_t * const p0 = &v1->p;

		// if v2 already exists in the edges, discard it
		if (v1->edges[i] == v2)
			return;
	}

	// if we reach this point, we need to insert the edge
	v1->edges[v1->num_edges++] = v2;
}


stl_vertex_t *
stl_vertex_find(
	stl_vertex_t ** const vertices,
	int * num_vertex_ptr,
	const v3_t * const p
)
{
	const int num_vertex = *num_vertex_ptr;

	for (int x = 0 ; x < num_vertex ; x++)
	{
		stl_vertex_t * const v = vertices[x];

		if (v3_eq(&v->p, p))
			return v;
	}

	if (debug)
	fprintf(stderr, "%d: %f,%f,%f\n",
		num_vertex,
		p->p[0],
		p->p[1],
		p->p[2]
	);

	stl_vertex_t * const v = vertices[(*num_vertex_ptr)++] = calloc(1, sizeof(*v));
	v->p = *p;
	return v;
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
	const float thick = 1;
	const int do_square = 0;

	fprintf(stderr, "header: '%s'\n", hdr->header);
	fprintf(stderr, "num: %d\n", num_triangles);

	// generate the unique list of vertices and their
	// correponding edges
	stl_vertex_t ** const vertices = calloc(3*num_triangles, sizeof(*vertices));

	int num_vertex = 0;

	for(int i = 0 ; i < num_triangles ; i++)
	{
		stl_vertex_t * vp[3] = {};

		for (int j = 0 ; j < 3 ; j++)
		{
			const v3_t * const p = &stl_faces[i].p[j];
			vp[j] = stl_vertex_find(vertices, &num_vertex, p);
		}

		// walk all of other triangles to figure out if
		// any of the triangles are coplanar and have shared
		// edges.

		uint8_t coplanar_mask = 0;
		for (int j = 0 ; j < num_triangles ; j++)
		{
			if (j == i)
				continue;

			coplanar_mask |= coplanar_check(
				&stl_faces[i], &stl_faces[j]);
		}

		// all three vertices are mapped; generate the
		// connections
		for (int j = 0 ; j < 3 ; j++)
		{
			stl_vertex_t * const v = vp[j];

			// if the edge from j to j+1 is not coplanar,
			// add it to the list
			if ((coplanar_mask & (1 << j)) == 0)
				stl_edge_insert(v, vp[(j+1) % 3]);

			// if the edge from j+2 to j is not coplanar
			const uint8_t j2 = (j + 2) % 3;
			if ((coplanar_mask & (1 << j2)) == 0)
				stl_edge_insert(v, vp[j2]);
		}
	}

	fprintf(stderr, "%d unique vertices\n", num_vertex);
	for (int i = 0 ; i < num_vertex ; i++)
	{
		stl_vertex_t * const v = vertices[i];
		printf("translate([%f,%f,%f]) {\n",
			v->p.p[0],
			v->p.p[1],
			v->p.p[2]
		);
		
		printf("sphere(r=4); // %p\n", v);

		for (int j = 0 ; j < v->num_edges ; j++)
		{
			stl_vertex_t * const v2 = v->edges[j];
			const v3_t d = v3_sub(v2->p, v->p);
			const float len = v3_len(&v2->p, &v->p);

			const float b = acos(d.p[2] / len) * 180/M_PI;
			const float c = d.p[0] == 0 ? sign(d.p[1]) * 90 : atan2(d.p[1], d.p[0]) * 180/M_PI;
//
			printf("%%rotate([0,%f,%f]) cylinder(r=1, h=%f); // %p\n",
				b,
				c,
				len*.45,
				v2
			);
		}

		printf("}\n");
	}

	return 0;
}
