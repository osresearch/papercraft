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


typedef struct face face_t;
typedef struct poly poly_t;

struct face
{
	float sides[3];
	face_t * next[3];
	int next_edge[3];
	int coplanar[3];
	int used[3];
};

// once this triangle has been used, it will be placed
// in a polygon group and fixed in a position relative to that group
struct poly
{
	int start_edge;
	int printed;

	// local coordinates of the triangle vertices
	float a;
	float x2;
	float y2;
	float rot;

	// absolute coordintes of the triangle vertices
	float p[3][2];

	// todo: make this const and add backtracking
	face_t * face;
	poly_t * next[3];

	poly_t * work_next;
};


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



/* Returns the 0 for coplanar, negative for mountain, positive for valley.
 * (approximates the angle between two triangles that share one edge).
 */
int
coplanar_check(
	const stl_face_t * const f1,
	const stl_face_t * const f2
)
{
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
	if (debug) fprintf( stderr, "%p %p %s: %f\n", f1, f2, check ? "yes" : "no", dot);
	return (int) dot;
}


/** Translate a list of STL triangles into a connected graph of faces.
 *
 * If there are any triangles that do not have three connected edges,
 * the first error will be reported and NULL will be returned.
 */
face_t *
stl2faces(
	const stl_face_t * const stl_faces,
	const int num_triangles
)
{
	face_t * const faces = calloc(num_triangles, sizeof(*faces));

	// convert the stl triangles into faces
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];
		face_t * const f = &faces[i];

		f->sides[0] = v3_len(&stl->p[0], &stl->p[1]);
		f->sides[1] = v3_len(&stl->p[1], &stl->p[2]);
		f->sides[2] = v3_len(&stl->p[2], &stl->p[0]);
		if (debug) fprintf(stderr, "%p %f %f %f\n",
			f, f->sides[0], f->sides[1], f->sides[2]);
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
		fprintf(stderr, "%d missing edges?\n", i);
		free(faces);
		return NULL;
	}

	return faces;
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

		// all three vertices are mapped; generate the
		// connections
		for (int j = 0 ; j < 3 ; j++)
		{
			stl_vertex_t * const v = vp[j];
			stl_vertex_find(v->edges, &v->num_edges, &vp[(j+1) % 3]->p);
			stl_vertex_find(v->edges, &v->num_edges, &vp[(j+2) % 3]->p);
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
				len/3,
				v2
			);
		}

		printf("}\n");
	}

#if 0

	face_t * const faces = stl2faces(stl_faces, num_triangles);

	for (int i = 0 ; i < num_triangles ; i++)
	{
		//if (i > 20) break;
		const stl_face_t * const raw = &stl_faces[i];
		face_t * const f = &faces[i];
		for (int j = 0 ; j < 3 ; j++)
		{
			// if this edge is coplanar with the other
			// triangle, do not output it.
			if (f->coplanar[j] == 0)
				continue;

			// if we have already transited this edge
			if (f->used[j])
				continue;

			// flag that we have used this vertex on the adjacent

			f->used[j] = 1;

			const v3_t * const p1 = &raw->p[(j+0) % 3];
			const v3_t * const p2 = &raw->p[(j+1) % 3];
			const float len = v3_len(p1,p2);
			const v3_t d = v3_sub(*p2, *p1);
			if (len == 0)
				continue;

			printf("translate([%f,%f,%f]) sphere(r=%f);\n",
				p1->p[0], p1->p[1], p1->p[2],
				4*thick
			);
				

			const float b = acos(d.p[2] / len) * 180/M_PI;
			const float c = d.p[0] == 0 ? sign(d.p[1]) * 90 : atan2(d.p[1], d.p[0]) * 180/M_PI;
//
			// generate a cube that goes from
			// p1 to p2.
			printf("%%translate([%f,%f,%f]) rotate([%f,%f,%f]) ",
				p1->p[0], p1->p[1], p1->p[2],
				0.0, b, c
			);

			if (do_square)
			{
				printf("translate([0,0,%f]) cube([%f,%f,%f], center=true);\n",
					len/2,
					thick,
					thick,
					len
				);
			} else {
				printf("cylinder(r=%f, h=%f);\n",
					thick,
					len
				);
			}
		}
	}
#endif


	return 0;
}
