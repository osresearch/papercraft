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

	for(int i = 0 ; i < num_triangles ; i++)
	{
		const stl_3d_file_triangle_t * const ft = &fts[i];
		stl_face_t * const f = &stl->face[i];

		for (int j = 0 ; j < 3 ; j++)
		{
			const v3_t * const p = &ft->p[j];

			stl_vertex_t * const v
				= stl_vertex_find(stl->vertex, &stl->num_vertex, p);

			// add this vertex to this face
			f->vertex[j] = v;

			// and add this face to the vertex
			v->face[v->num_face++] = f;
		}
	}

	return stl;
}
