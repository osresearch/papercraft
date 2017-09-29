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
#include "v3.h"
#include "camera.h"

#ifndef M_PI
#define 	M_PI   3.1415926535897932384
#endif

static int debug = 1;

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
	v3_t p[3];
	v3_t normal;
	float area;
	float min[3];
	float max[3];
	tri_t * next;
	tri_t ** prev;
};


// line segment has to track its source so that it knows which to not
// compare against in its occlusion checks.
typedef struct _seg_t seg_t;
struct _seg_t {
	v3_t p[2];
	v3_t src[2];
	seg_t * next;
};


#if 0
typedef struct face face_t;
typedef struct poly poly_t;

struct face
{
	float sides[3];
	face_t * next[3];
	int next_edge[3];
	int coplanar[3];
	int used;
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
#endif


void
svg_line(
	const char * color,
	const float * p1,
	const float * p2,
	int dash
)
{
	if (!dash)
	{
		printf("<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%s\" stroke-width=\"0.5px\"/>\n",
			p1[0],
			p1[1],
			p2[0],
			p2[1],
			color
		);
		return;
	}

	// dashed line, split in the middle
	const float dx = p2[0] - p1[0];
	const float dy = p2[1] - p1[1];

	const float h1[] = {
		p1[0] + dx*0.45,
		p1[1] + dy*0.45,
	};
	const float h2[] = {
		p1[0] + dx*0.55,
		p1[1] + dy*0.55,
	};


	svg_line(color, p1, h1, 0);
	svg_line(color, h2, p2, 0);
}


#if 0
void
rotate(
	float * p,
	const float * origin,
	float a,
	float x,
	float y
)
{
	p[0] = cos(a) * x - sin(a) * y + origin[0];
	p[1] = sin(a) * x + cos(a) * y + origin[1];
}


/* Rotate and translate a triangle */
void
poly_position(
	poly_t * const g,
	const poly_t * const g_src,
	float rot,
	float trans_x,
	float trans_y
)
{
	const face_t * const f = g->face;
	const int start_edge = g->start_edge;

	float a = f->sides[(start_edge + 0) % 3];
	float c = f->sides[(start_edge + 1) % 3];
	float b = f->sides[(start_edge + 2) % 3];
	float x2 = (a*a + b*b - c*c) / (2*a);
	float y2 = sqrt(b*b - x2*x2);

	// translate by trans_x/trans_y in the original ref frame
	// to get the origin point
	float origin[2];
	rotate(origin, g_src->p[0], g_src->rot, trans_x, trans_y);

	g->rot = g_src->rot + rot;
	g->a = a;
	g->x2 = x2;
	g->y2 = y2;

//fprintf(stderr, "%p %d %f %f %f %f => %f %f %f\n", f, start_edge, g->rot*180/M_PI, a, b, c, x2, y2, rot);
	rotate(g->p[0], origin, g->rot, 0, 0);
	rotate(g->p[1], origin, g->rot, a, 0);
	rotate(g->p[2], origin, g->rot, x2, y2);
}


static void
enqueue(
	poly_t * g,
	poly_t * const new_g,
	int at_head
)
{
	if (at_head)
	{
		new_g->work_next = g->work_next;
		g->work_next = new_g;
		return;
	}

	// go to the end of the line
	while (g->work_next)
		g = g->work_next;
	g->work_next = new_g;
}


static poly_t * poly_root;
static float poly_min[2], poly_max[2];
#endif

static inline int
v2_eq(
	const float p0[],
	const float p1[]
)
{
	const float dx = p0[0] - p1[0];
	const float dy = p0[1] - p1[1];

	// are the points within epsilon of each other?
	if (-EPS < dx && dx < EPS
	&&  -EPS < dy && dy < EPS)
		return 1;

	// nope, not equal
	return 0;
}

static inline int
v2_dist(
	const float p0[],
	const float p1[]
)
{
	const float dx = p0[0] - p1[0];
	const float dy = p0[1] - p1[1];

	return sqrt(dx*dx + dy*dy);
}


// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
int
get_line_intersection(
	float p0_x,
	float p0_y,
	float p1_x,
	float p1_y, 
	float p2_x,
	float p2_y,
	float p3_x,
	float p3_y,
	float *i_x,
	float *i_y
)
{
	float s1_x = p1_x - p0_x;
	float s1_y = p1_y - p0_y;
	float s2_x = p3_x - p2_x;
	float s2_y = p3_y - p2_y;

	float s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y))
		/ (-s2_x * s1_y + s1_x * s2_y);

	float t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x))
		/ (-s2_x * s1_y + s1_x * s2_y);

	if (s > EPS && s < 1-EPS && t > EPS && t < 1-EPS)
	{
		if(0) fprintf(stderr, "collision: %f,%f->%f,%f %f,%f->%f,%f == %f,%f\n",
			p0_x, p0_y,
			p1_x, p1_y,
			p2_x, p2_y,
			p3_x, p3_y,
			s,
			t
		);

		// Collision detected
		if (i_x != NULL)
			*i_x = p0_x + (t * s1_x);
		if (i_y != NULL)
			*i_y = p0_y + (t * s1_y);
		return 1;
	}

	return 0; // No collision
}


int
intersect(
	const v3_t * const p00,
	const v3_t * const p01,
	const v3_t * const p10,
	const v3_t * const p11,
	float *px,
	float *py
)
{
	// special case; if this is the same line, it does not intersect
	if (v2_eq(p00->p, p10->p) && v2_eq(p01->p, p11->p))
		return 0;
	if (v2_eq(p01->p, p10->p) && v2_eq(p00->p, p11->p))
		return 0;

	return get_line_intersection(
		p00->p[0], p00->p[1],
		p01->p[0], p01->p[1],
		p10->p[0], p10->p[1],
		p11->p[0], p11->p[1],
		px,
		py
	);
}


#if 0
/** Check to see if two triangles overlap */
int
overlap_poly(
	const poly_t * const g1,
	const poly_t * const g2
)
{
	if (intersect(g1->p[0], g1->p[1], g2->p[0], g2->p[1]))
		return 1;
	if (intersect(g1->p[0], g1->p[1], g2->p[1], g2->p[2]))
		return 1;
	if (intersect(g1->p[0], g1->p[1], g2->p[2], g2->p[0]))
		return 1;

	if (intersect(g1->p[1], g1->p[2], g2->p[0], g2->p[1]))
		return 1;
	if (intersect(g1->p[1], g1->p[2], g2->p[1], g2->p[2]))
		return 1;
	if (intersect(g1->p[1], g1->p[2], g2->p[2], g2->p[0]))
		return 1;

	if (intersect(g1->p[2], g1->p[0], g2->p[0], g2->p[1]))
		return 1;
	if (intersect(g1->p[2], g1->p[0], g2->p[1], g2->p[2]))
		return 1;
	if (intersect(g1->p[2], g1->p[0], g2->p[2], g2->p[0]))
		return 1;

	return 0;
}


/** Check to see if any triangles overlap */
int
overlap_check(
	const poly_t * g,
	const poly_t * const new_g
)
{
	// special case -- if the root is the same as the one that we
	// are checking, then it does not overlap
	if (g == new_g)
		return 0;

	while (g)
	{
		if (overlap_poly(g, new_g))
			return 1;
		g = g->work_next;
	}

	return 0;
}


/** recursively try to fix up the triangles.
 *
 * returns the maximum number of triangles added
 */
int
poly_build(
	poly_t * const g
)
{
	face_t * const f = g->face;
	const int start_edge = g->start_edge;
	f->used = 1;

	// update the group's bounding box
	for (int i = 0 ; i < 3 ; i++)
	{
		const float px = g->p[i][0];
		const float py = g->p[i][1];

		if (px < poly_min[0]) poly_min[0] = px;
		if (px > poly_max[0]) poly_max[0] = px;

		if (py < poly_min[1]) poly_min[1] = py;
		if (py > poly_max[1]) poly_max[1] = py;
	}
		

	if (debug) fprintf(stderr, "%p: adding to poly\n", f);

   for(int pass = 0 ; pass < 2 ; pass++)
   {
	// for each edge, find the triangle that matches
	for (int i = 0 ; i < 3 ; i++)
	{
		const int edge = (i + start_edge) % 3;
		face_t * const f2 = f->next[edge];
		assert(f2 != NULL);
		if (f2->used)
			continue;
		if (pass == 0 && f->coplanar[edge] == 0)
			continue;

		// create a group that translates and rotates
		// such that it lines up with this edge
		float trans_x, trans_y, rotate;
		if (i == 0)
		{
			trans_x = g->a;
			trans_y = 0;
			rotate = M_PI;
		} else
		if (i == 1)
		{
			trans_x = g->x2;
			trans_y = g->y2;
			rotate = -atan2(g->y2, g->a - g->x2);
		} else
		if (i == 2)
		{
			trans_x = 0;
			trans_y = 0;
			rotate = atan2(g->y2, g->x2);
		} else {
			errx(EXIT_FAILURE, "edge %d invalid?\n", i);
		}

		// position this one translated and rotated
		poly_t * const g2 = calloc(1, sizeof(*g2));
		g2->face = f2;
		g2->start_edge = f->next_edge[edge];

		poly_position(
			g2,
			g,
			rotate, 
			trans_x,
			trans_y
		);

		if (overlap_check(poly_root, g2))
		{
			free(g2);
			continue;
		}

		// no overlap, add it to the current group
		g->next[i] = g2;
		g2->next[0] = g;
		f2->used = 1;

		// if g2 is a coplanar triangle, process it now rather than
		// defering the work.
		if (f->coplanar[edge] == 0)
			enqueue(g, g2, 1);
		else
			enqueue(g, g2, 0);
	}
    }

	return 0;
}


void
svg_text(
	float x,
	float y,
	float angle,
	const char * fmt,
	...
)
{

	printf("<g transform=\"translate(%f %f) rotate(%f)\">",
		x,
		y,
		angle
	);

	printf("<text x=\"-2\" y=\"1.5\" style=\"font-size:1.5px;\">");

	va_list ap;
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);

	printf("</text></g>\n");
}

void
poly_print(
	poly_t * const g
)
{
	const face_t * const f = g->face;
	const int start_edge = g->start_edge;

	g->printed = 1;

	// draw this triangle;
	// if the edge is an outside, which means that the group
	// has no next element, draw a cut line.  If there is an
	// adjacent neighbor and it is not coplanar, draw a score line
printf("<g><!-- %p %d %f %f->%p %f->%p %f->%p -->\n",
	f,
	g->start_edge, g->rot * 180/M_PI,
	f->sides[0],
	f->next[0],
	f->sides[1],
	f->next[1],
	f->sides[2],
	f->next[2]
);

	int cut_lines = 0;
	const uintptr_t a1 = (0x7FFFF & (uintptr_t) f) >> 3;

	for (int i = 0 ; i < 3 ; i++)
	{
		const int edge = (start_edge + i) % 3;
		poly_t * const next = g->next[i];

		if (!next)
		{
			// draw a cut line
			const float * const p1 = g->p[i];
			const float * const p2 = g->p[(i+1) % 3];
			const float cx = (p2[0] + p1[0]) / 2;
			const float cy = (p2[1] + p1[1]) / 2;
			const float dx = (p2[0] - p1[0]);
			const float dy = (p2[1] - p1[1]);
			const float angle = atan2(dy, dx) * 180 / M_PI;

			svg_line("#FF0000", p1, p2, 0);
			cut_lines++;

			// use the lower address as the label
			if (draw_labels)
			{
				uintptr_t a2 = (0x7FFFF & (uintptr_t) f->next[edge]) >> 3;
				if (a2 > a1)
					a2 = a1;
				svg_text(cx, cy, angle, "%04x", a2);
			}

			continue;
		}

		if (next->printed)
			continue;

		if (f->coplanar[edge] < 0)
		{
			// draw a mountain score line since they are not coplanar
			svg_line("#00FF00", g->p[i], g->p[(i+1) % 3], 1);
		} else
		if (f->coplanar[edge] > 0)
		{
			// draw a valley score line since they are not coplanar
			svg_line("#00FF00", g->p[i], g->p[(i+1) % 3], 0);
		} else {
			// draw a shadow line since they are coplanar
			//svg_line("#F0F0F0", g->p[i], g->p[(i+1) % 3]);
		}
	}

/*
	// only draw labels if requested and if there are any cut-edges
	// on this polygon.
	const float tx = (g->p[0][0] + g->p[1][0] + g->p[2][0]) / 3.0;
	const float ty = (g->p[0][1] + g->p[1][1] + g->p[2][1]) / 3.0;
	if (draw_labels && cut_lines > 0)
	svg_text(tx, ty, 0, "%04x",
		(0x7FFFF & (uintptr_t) f) >> 3);
*/

printf("</g>\n");

	for (int i = 0 ; i < 3 ; i++)
	{
		poly_t * const next = g->next[i];
		if (!next || next->printed)
			continue;

		poly_print(next);
	}
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
#endif

/*
s = 1/(2*Area)*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
t = 1/(2*Area)*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);
where Area is the (signed) area of the triangle:

Area = 0.5 *(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
Just evaluate s, t and 1-s-t. The point p is inside the triangle if and only if they are all positive.
*/
int
tri_inside(
	const tri_t * const t,
	const v3_t * const p
)
{
	const float p0x = t->p[0].p[0];
	const float p0y = t->p[0].p[1];
	const float p1x = t->p[1].p[0];
	const float p1y = t->p[1].p[1];
	const float p2x = t->p[2].p[0];
	const float p2y = t->p[2].p[1];

	const float px = p->p[0];
	const float py = p->p[1];

	const float u = p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py;
	const float v = p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py;

	if (u <= 0 || v <= 0)
		return 0;

	// maybe inside; check for sure
	if (u + v >= 2 * t->area)
		return 0;

	// inside!
if(0) fprintf(stderr, "%p: %f,%f inside %f,%f %f,%f %f,%f\n",
	t,
	px, py,
	p0x, p0y,
	p1x, p1y,
	p2x, p2y
);

	return 1;
}


tri_t *
tri_new(
	const v3_t * p
)
{
	tri_t * const t = calloc(1, sizeof(*t));
	if (!t)
		return NULL;
	for(int i = 0 ; i < 3  ; i++)
		t->p[i] = p[i];

	// precompute the area
	const float p0x = t->p[0].p[0];
	const float p0y = t->p[0].p[1];
	const float p1x = t->p[1].p[0];
	const float p1y = t->p[1].p[1];
	const float p2x = t->p[2].p[0];
	const float p2y = t->p[2].p[1];

	t->area = 0.5 *(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);

	// precompute the normal
	t->normal = v3_cross(
		v3_sub(t->p[1], t->p[0]),
		v3_sub(t->p[2], t->p[1])
	);


	// compute the bounding box for the triangle
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


int
tri_occluded(
	const tri_t * zlist,
	const tri_t * t
)
{
	for( const tri_t * t2 = zlist ; t2 ; t2 = t2->next )
	{
		if (t2 == t)
			continue;

		// if any of the points of t are outside of t2,
		// then t2 does not totally occlude t
		if (!tri_inside(t2, &t->p[0]))
			continue;
		if (!tri_inside(t2, &t->p[1]))
			continue;
		if (!tri_inside(t2, &t->p[2]))
			continue;

		// if any point of t2 is behind t, then it does not occlude
		// (might intersect, but we don't handle that)
		if (t2->min[2] > t->min[2])
			continue;

		// looks like we might be occluded
		return 1;
	}

	// probably not occluded
	return 0;
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


/*
int
tri_line_intersect(
	const tri_t * zlist,
	const tri_t * t
)
{
	for( const tri_t * t2 = zlist ; t2 ; t2 = t2->next )
	{
		if (t2 == t)
			continue;
		for(int j = 0 ; j < 3 ; j++)
		{
			const v3_t * const p0 = &t->p[j].p;
			const v3_t * const p1 = &t->p[(j+1) % 3].p;
*/


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
	const float p0x = s->p[0].p[0];
	const float p0y = s->p[0].p[1];
	const float p0z = s->p[0].p[2];
	const float p1x = s->p[1].p[0];
	const float p1y = s->p[1].p[1];
	const float p1z = s->p[1].p[2];

	for( const tri_t * t = zlist ; t ; t = t->next )
	{
		// if the segment is closer than the triangle,
		// then we no longer have to check any further into
		// the zlist (it is sorted by depth).
		if (p0z < t->min[2] && p1z < t->min[2])
			break;

#if 0
		// do a quick test of does this segment even comes
		// close to this triangle
		if (p0x < t->min[0] && p1x < t->min[0]
		&&  p0y < t->min[1] && p1y < t->min[1])
			continue;
		if (p0x > t->max[0] && p1x > t->max[0]
		&&  p0y > t->max[2] && p1y > t->max[2])
			continue;
		if (p0x < t->min[0] && p1x < t->min[0]
		&&  p0y > t->max[2] && p1y > t->max[2])
			continue;
		if (p0x > t->max[0] && p1x > t->max[0]
		&&  p0y < t->min[1] && p1y < t->min[1])
			continue;
		
		// make sure this isn't the same actual line
		if (v3_eq(&s->src[0], &t->p[0]) || v3_eq(&s->src[1], &t->p[1]))
			continue;
		if (v3_eq(&s->src[0], &t->p[1]) || v3_eq(&s->src[1], &t->p[2]))
			continue;
		if (v3_eq(&s->src[0], &t->p[2]) || v3_eq(&s->src[1], &t->p[0]))
			continue;
#endif

		int inside0 = tri_inside(t, &s->p[0]);
		int inside1 = tri_inside(t, &s->p[1]);

		// if both are inside we discard this segment
		if (inside0 && inside1)
			return;

		// split the segment for each intersection with the
		// triangle segments and add it to the work queue.
		int intersections = 0;
		v3_t ix[3] = {};
 		const float max_z = max(s->p[0].p[2], s->p[1].p[2]);

		for(int j = 0 ; j < 3 ; j++)
		{
			ix[j].p[2] = max_z;
			int rc = intersect(
				&s->p[0], &s->p[1],
				&t->p[j], &t->p[(j+1)%3],
				&ix[intersections].p[0], &ix[intersections].p[1]
			);

			if (!rc)
				continue;

			intersections++;
		}

		// if none of them intersect, we keep looking
		if (intersections == 0)
			continue;

fprintf(stderr, "split %d %d inter %d\n", inside0 , inside1, intersections);
		if (intersections == 3)
		{
			fprintf(stderr, "uh, three intersections?\n");
			return;
		}

		if (intersections == 2)
		{
			if (inside0 || inside1)
			{
				fprintf(stderr, "uh, inside but two intersections?\n");
				return;
			}

			// we have to create a new segment
			// and shorten the existing segment
			// find the two intersections that we have
			// update the src field
			

			fprintf(stderr, "two intersections\n");
			const float d0 = v2_dist(s->p[0].p, ix[0].p);
			const float d1 = v2_dist(s->p[1].p, ix[0].p);
			seg_t * news;
			if (d0 < d1)
			{
				// split from p0 to ix0
				news = seg_new(s->p[0], ix[0]);
				news->src[1] = s->p[1];
				s->p[1] = ix[1];
			} else {
				// split from p0 to ix1
				news = seg_new(s->p[0], ix[1]);
				news->src[1] = s->p[1];
				s->p[1] = ix[0];
			}

			// recursively start splitting the new segment
			// starting at our current z-depth
			tri_seg_intersect(zlist, news, slist_visible);

			// continue splitting our current segment
			continue;
		}

		if (intersections == 1)
		{
fprintf(stderr, "split %d %d\n", inside0, inside1);
			if (inside0)
			{
				// shorten it on the 0 side
				s->p[0] = ix[0];
			} else
			if (inside1)
			{
				// shorten it on the 1 side
				s->p[1] = ix[0];
			} else {
				fprintf(stderr, "uh, both outside but one intersection?\n");
				return;
			}
		}
			

if(0) fprintf(stderr, "check: %.0f,%.0f -> %.0f,%.0f  %.0f,%.0f %.0f,%.0f %.0f,%.0f\n",
	s->p[0].p[0],
	s->p[0].p[1],
	s->p[1].p[0],
	s->p[1].p[1],
	t->p[0].p[0],
	t->p[0].p[1],
	t->p[1].p[0],
	t->p[1].p[1],
	t->p[2].p[0],
	t->p[2].p[1]
);
		//return;
	}

	// if we've reached here the segment is visible
	// and should be added to the visible list
if(0) fprintf(stderr, "good: %.0f,%.0f,%.0f-> %.0f,%.0f,%.0f\n",
	s->p[0].p[0],
	s->p[0].p[1],
	s->p[0].p[2],
	s->p[1].p[0],
	s->p[1].p[1],
	s->p[1].p[2]
);

	s->next = *slist_visible;
	*slist_visible = s;
}

			

int main(
	int argc,
	char ** argv
)
{
	const size_t max_len = 32 << 20;
	uint8_t * const buf = calloc(max_len, 1);

	float phi = argc > 1 ? atof(argv[1]) * M_PI/180 : 0;
	float theta = argc > 2 ? atof(argv[2]) * M_PI/180 : 0;
	float psi = argc > 3 ? atof(argv[3]) * M_PI/180 : 0;

	ssize_t rc = read(0, buf, max_len);
	if (rc == -1)
		return EXIT_FAILURE;

	const stl_header_t * const hdr = (const void*) buf;
	const stl_face_t * const stl_faces = (const void*)(hdr+1);
	const int num_triangles = hdr->num_triangles;

	if(debug)
	{
		fprintf(stderr, "header: '%s'\n", hdr->header);
		fprintf(stderr, "num: %d\n", num_triangles);
	}


	// looking at (0,0,0)
	v3_t eye = { { 0, 0, 400 } };
	const camera_t * const cam = camera_new(eye, phi, theta, psi);

	printf("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");

	float off_x = 500;
	float off_y = 500;
	printf("<g transform=\"translate(%f %f)\">\n", off_x, off_y);

	int rejected = 0;
	tri_t * zlist = NULL;
	seg_t * slist = NULL;
	seg_t * slist_visible = NULL;

	// transform the stl by the camera projection and generate
	// a z-sorted list of triangles
	for (int i = 0 ; i < num_triangles ; i++)
	{
		const stl_face_t * const stl = &stl_faces[i];

		v3_t s[3];

		for(int j = 0 ; j < 3 ; j++)
			camera_project(cam, &stl->p[j], &s[j]);

		tri_t * const tri = tri_new(s);

		// reject this face if any of the vertices are behind us
		if (tri->min[2] < 0)
			goto reject;

		// do a back-face cull to determine if this triangle
		// is not facing us. we have to determine the orientation
		// from the winding of the new projection
		if (tri->normal.p[2] <= 0)
			goto reject;

		// it passes the first tests, so insert the triangle
		// into the list and the three line segments
		tri_insert(&zlist, tri);

		for(int j = 0 ; j < 3 ; j++)
		{
			seg_t * s = seg_new(tri->p[j], tri->p[(j+1) % 3]);
			s->next = slist;
			slist = s;
		}

		continue;

reject:
		tri_delete(tri);
		rejected++;
	}

	if (debug)
		fprintf(stderr, "Rejected %d triangles\n", rejected);

	// we now have a z-sorted list of triangles
	rejected = 0;

	// work on each segment, intersecting it with all of the triangles
	while(slist)
	{
		seg_t * s = slist;
		slist = s->next;

		tri_seg_intersect(zlist, s, &slist_visible);

	}

	// display all of the visible segments
	for(seg_t * s = slist_visible ; s ; s = s->next)
	{
		svg_line("#FF0000", s->p[0].p, s->p[1].p, 0);
	}

	if (debug)
		fprintf(stderr, "Occluded %d triangles\n", rejected);


	printf("</g>\n");
	printf("</svg>\n");

	return 0;
}
