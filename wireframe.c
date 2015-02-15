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


#define MAX_VERTEX 64

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

	if (debug)
		fprintf(stderr, "%p %p: %d %d %d\n",
			f1,
			f2,
			match[0], match[1], match[2]
		);

	// otherwise they do not share enough points
	if (mask == 0)
		return 0;

#if 0
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

	if (debug)
	fprintf(stderr, "dot %f:\n %f,%f,%f\n %f,%f,%f\n %f,%f,%f\n %f,%f,%f\n",
		dot,
		x1.p[0], x1.p[1], x1.p[2],
		x2.p[0], x2.p[1], x2.p[2],
		x3.p[0], x3.p[1], x3.p[2],
		x4.p[0], x4.p[1], x4.p[2]
	);
	
	int check = -EPS < dot && dot < +EPS;

	// if the dot product is not close enough to zero, they
	// are not coplanar.
	if (!check)
		return 0;

	// coplanar! return the shared edge mask
	return mask;
#else
	// if the normals are close enough, then it is coplanner
	if (v3_eq(&f1->normal, &f2->normal))
		return mask;
	else
		return 0;
#endif
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
		// if v2 already exists in the edges, discard it
		if (v1->edges[i] == v2)
			return;
	}

	// if we reach this point, we need to insert the edge
	if (debug)
	fprintf(stderr, "%p: edge %d -> %p\n",
		v1,
		v1->num_edges,
		v2
	);
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

	if (num_vertex == 6)
	fprintf(stderr, "%f %f %f\n",
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
	const float thick = 7.8;
	const int do_square = 1;

	fprintf(stderr, "header: '%s'\n", hdr->header);
	fprintf(stderr, "num: %d\n", num_triangles);

	// generate the unique list of vertices and their
	// correponding edges
	stl_vertex_t ** const vertices = calloc(3*num_triangles, sizeof(*vertices));

	int num_vertex = 0;

	for(int i = 0 ; i < num_triangles ; i++)
	{
		if (debug) fprintf(stderr, "---------- triangle %d (%d)\n", i, num_vertex);

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

			if (debug)
			fprintf(stderr, "%d: check %d -> %d\n", num_vertex, i, j);

			coplanar_mask |= coplanar_check(
				&stl_faces[i], &stl_faces[j]);
		}

		if (debug)
			fprintf(stderr, "mask %d\n", coplanar_mask);

		// all three vertices are mapped; generate the
		// connections
		for (int j = 0 ; j < 3 ; j++)
		{
			stl_vertex_t * const v = vp[j];

			// if the edge from j to j+1 is not coplanar,
			// add it to the list
			if ((coplanar_mask & (1 << j)) == 0)
			{
				if (debug)
				fprintf(stderr, "%p: %d insert\n", v, j);
				stl_edge_insert(v, vp[(j+1) % 3]);
			}

/*
			// if the edge from j+2 to j is not coplanar
			const uint8_t j2 = (j + 2) % 3;
			if ((coplanar_mask & (1 << j2)) == 0)
			{
				if (debug)
				fprintf(stderr, "%p: %d insert back\n", v, j2);
				stl_edge_insert(v, vp[j2]);
			}
*/
		}
	}

	fprintf(stderr, "%d unique vertices\n", num_vertex);
	printf("thick=%f;\n"
		"module connector(len) {\n"
		"  render() difference() {\n"
		"    cylinder(r=thick/2+2, h=2*thick);\n"
		//"    translate([0,0,len/2+2]) cube([thick,thick,2*thick]);\n"
		"    translate([0,0,thick/2+2]) cylinder(r=thick/2, h=2*thick);\n"
		"  }\n"
		//"  %%translate([0,0,len*0.48/2]) cube([thick,thick,len*0.48], center=true);\n"
		"  %%translate([0,0,0]) cylinder(r=thick/2, h=len*0.48);\n"
		"}\n",
		thick
	);

	for (int i = 0 ; i < num_vertex ; i++)
	{
		stl_vertex_t * const v = vertices[i];
		printf("translate([%f,%f,%f]) {\n",
			v->p.p[0],
			v->p.p[1],
			v->p.p[2]
		);
		
		printf("sphere(r=%f); // %d %p\n", thick/2+2, i, v);

		for (int j = 0 ; j < v->num_edges ; j++)
		{
			stl_vertex_t * const v2 = v->edges[j];
			const v3_t d = v3_sub(v2->p, v->p);
			const float len = v3_len(&v2->p, &v->p);

			const float b = acos(d.p[2] / len) * 180/M_PI;
			const float c = d.p[0] == 0 ? sign(d.p[1]) * 90 : atan2(d.p[1], d.p[0]) * 180/M_PI;
//
			printf("rotate([0,%f,%f]) ", b, c);

			if (do_square)
				printf("connector(%f);\n", len);
			else
				printf(" cylinder(r=1, h=%f); // %p\n",
					len*.45,
					v2
				);
		}

		printf("}\n");
	}

	return 0;
}
