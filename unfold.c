#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>


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
	v3_t p0;
	v3_t p1;
	v3_t p2;
} __attribute__((__packed__))
stl_face_t;


#define MAX_POINTS 24

typedef struct
{
	int n;
	int p[MAX_POINTS];
} poly_t;


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


int main(void)
{
	const size_t max_len = 1 << 20;
	uint8_t * const buf = calloc(max_len, 1);

	ssize_t rc = read(0, buf, max_len);
	if (rc == -1)
		return EXIT_FAILURE;

	const stl_header_t * const hdr = (const void*) buf;
	const stl_face_t * const faces = (const void*)(hdr+1);
	const int num_triangles = hdr->num_triangles;

	fprintf(stderr, "header: '%s'\n", hdr->header);
	fprintf(stderr, "num: %d\n", num_triangles);

	// worst case -- all separate polygons
	poly_t * const polys = calloc(num_triangles, sizeof(*polys));

	// Collapse coplanar triangles into larger polygons
	v3_t * const vertices = calloc(num_triangles, 3);
	int num_vertices = 0;

	for(int i = 0 ; i < num_triangles ; i++)
	{
		// see if this matches an existing vertex
		const stl_face_t * const t = &faces[i];
		poly_t * const p = &polys[i];
		p->n = 3;
		p->p[0] = p->p[1] = p->p[2] = -1;

		for (int j = 0 ; j < num_vertices ; j++)
		{
			const v3_t * const v = &vertices[j];
			if (p->p[0] == -1 && v3_eq(v, &t->p0))
				p->p[1] = j;
			if (p->p[1] == -1 && v3_eq(v, &t->p1))
				p->p[1] = j;
			if (p->p[2] == -1 && v3_eq(v, &t->p2))
				p->p[2] = j;

			// check if we've found all of them
			if (p->p[0] >= 0 && p->p[1] >= 0 && p->p[2] >= 0)
				break;
		}

		// create new points if we haven't found matches
		if (p->p[0] < 0)
		{
			p->p[0] = num_vertices;
			vertices[num_vertices++] = t->p0;
		}
		if (p->p[1] < 0)
		{
			p->p[1] = num_vertices;
			vertices[num_vertices++] = t->p1;
		}
		if (p->p[3] < 0)
		{
			p->p[3] = num_vertices;
			vertices[num_vertices++] = t->p2;
		}
	}

	fprintf(stderr, "unique vertices: %d\n", num_vertices);

	return 0;
}
