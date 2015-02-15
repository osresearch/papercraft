#include "stl_3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>


static const int debug = 0;

typedef struct
{
	char header[80];
	uint32_t num_triangles;
} __attribute__((__packed__))
stl_3d_file_header_t;


typedef struct
{
	v3_t normal;
	v3_t p[3];
	uint16_t attr;
} __attribute__((__packed__))
stl_3d_file_triangle_t;


/** Find or create a vertex */
static stl_vertex_t *
stl_vertex_find(
	stl_vertex_t * const vertices,
	int * num_vertex_ptr,
	const v3_t * const p
)
{
	const int num_vertex = *num_vertex_ptr;

	for (int x = 0 ; x < num_vertex ; x++)
	{
		stl_vertex_t * const v = &vertices[x];

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

	stl_vertex_t * const v = &vertices[(*num_vertex_ptr)++];
	v->p = *p;

	return v;
}


/** Check to see if the two faces share an edge.
 * \return 0 if no common edge, 1 if there is a shared link
 */
static int
stl_has_edge(
	const stl_face_t * const f,
	const stl_vertex_t * const v1,
	const stl_vertex_t * const v2
)
{
	if (f->vertex[0] != v1 && f->vertex[1] != v1 && f->vertex[2] != v1)
		return 0;
	if (f->vertex[0] != v2 && f->vertex[1] != v2 && f->vertex[2] != v2)
		return 0;

	return 1;
}


/** Compute the angle between the two planes.
 * This is an approximation:
 * \return 0 == coplanar, negative == valley, positive == mountain.
 */
static double
stl_angle(
	const stl_face_t * const f1,
	const stl_face_t * const f2
)
{
	// find the four distinct points
	v3_t x1 = f1->vertex[0]->p;
	v3_t x2 = f1->vertex[1]->p;
	v3_t x3 = f1->vertex[2]->p;
	v3_t x4;

	for (int i = 0 ; i < 3 ; i++)
	{
		x4 = f2->vertex[i]->p;
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
	if (check)
		return 0;

	if (dot < 0)
		return -1;
	else
		return +1;
}


static void
stl_find_neighbors(
	stl_3d_t * const stl,
	stl_face_t * const f1
)
{
	for(int i = 0 ; i < 3 ; i++)
	{
		const stl_vertex_t * const v1 = f1->vertex[(i+0) % 3];
		const stl_vertex_t * const v2 = f1->vertex[(i+1) % 3];

		for(int j = 0 ; j < stl->num_face ; j++)
		{
			stl_face_t * const f2 = &stl->face[j];

			// skip this triangle against itself
			if (f1 == f2)
				continue;

			// find if these two triangles share an edge
			if (!stl_has_edge(f2, v1, v2))
				continue;

			f1->face[i] = f2;
			f1->angle[i] = stl_angle(f1, f2);
		}
	}
}


stl_3d_t *
stl_3d_parse(
	int fd
)
{
	ssize_t rc;
	stl_3d_file_header_t hdr;

	rc = read(fd, &hdr, sizeof(hdr));
	if (rc != sizeof(hdr))
		return NULL;

	const int num_triangles = hdr.num_triangles;
	fprintf(stderr, "%d triangles\n", num_triangles);

	stl_3d_file_triangle_t * fts;
	const size_t file_len = num_triangles * sizeof(*fts);
 	fts = calloc(1, file_len);

	rc = read(fd, fts, file_len);
	if (rc < 0 || (size_t) rc != file_len)
		return NULL;

	stl_3d_t * const stl = calloc(1, sizeof(*stl));

	*stl = (stl_3d_t) {
		.num_vertex = 0,
		.num_face = num_triangles,
		.vertex = calloc(num_triangles, sizeof(*stl->vertex)),
		.face = calloc(num_triangles, sizeof(*stl->face)),
	};

	// build the unique set of vertices and their connection
	// to each face.
	for(int i = 0 ; i < num_triangles ; i++)
	{
		const stl_3d_file_triangle_t * const ft = &fts[i];
		stl_face_t * const f = &stl->face[i];

		for (int j = 0 ; j < 3 ; j++)
		{
			const v3_t * const p = &ft->p[j];

			stl_vertex_t * const v = stl_vertex_find(
				stl->vertex,
				&stl->num_vertex,
				p
			);

			// add this vertex to this face
			f->vertex[j] = v;

			// and add this face to the vertex
			v->face[v->num_face] = f;
			v->face_num[v->num_face] = j;
			v->num_face++;
		}
	}

	// build the connections between each face
	for(int i = 0 ; i < num_triangles ; i++)
	{
		stl_face_t * const f = &stl->face[i];
		stl_find_neighbors(stl, f);
	}

	return stl;
}
