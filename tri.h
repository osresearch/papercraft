/*
 * Triangle manipulations
 */
#ifndef _tri_h_
#define _tri_h_

#include "v3.h"
#include "seg.h"
#include "svg.h"

extern int tri_debug;

typedef struct _tri_t tri_t;
struct _tri_t
{
	v3_t p[3]; // camera space
	v3_t normal; // camera space
	v3_t normal_xyz; // original xyz space
	v3_t min; // camera space
	v3_t max; // camera space
	tri_t * next;
	tri_t ** prev;
};


tri_t *
tri_new(
	const v3_t * p_cam,
	const v3_t * p_xyz
);


// insert a triangle into our z-sorted list
void
tri_insert(
	tri_t ** zlist,
	tri_t * t
);


void
tri_delete(tri_t * t);


// Compute the 2D area of a triangle in screen space
// using Heron's formula
float
tri_area_2d(
	const tri_t * const t
);


void
tri_print(
	const tri_t * const t
);


/* Check if two triangles are coplanar and share an edge.
 *
 * Returns -1 if not coplanar, 0-2 for the edge in t0 that they share.
 */
int
tri_coplanar(
	const tri_t * const t0,
	const tri_t * const t1,
	const float coplanar_eps
);


/*
 * Find the Z point of an XY coordinate in a triangle.
 *
 * p can be written as a combination of t01 and t02,
 * p - t0 = a * (t1 - t0) + b * (t2 - t0)
 * setting t0 to 0, this becomes:
 * p = a * t1 + b * t2
 * which is two equations with two unknowns
 *
 * Returns true if the point is inside the triangle
 */
int
tri_find_z(
	const tri_t * const t,
	const v3_t * const p,
	float * const zout
);


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
);


/*
 * Given a line segment and a list of triangles,
 * find if the line segment crosses any triangle.
 * If it crosses a triangle the segment will be shortened
 * and an additional one might be created.
 * Recusively try intersecting the new segment (starting at the same triangle)
 * and then continue trying the shortened segment.
 *
 * Line segments will be added to the visible list.
 * Returns the number of new elements created
 */
int
tri_seg_hidden(
	const tri_t * zlist,
	seg_t * s,
	seg_t ** slist_visible
);


/*
 * Fast check to see if t2 is entire occluded by t.
 */
int
tri_behind(
	const tri_t * const t,
	const tri_t * const t2
);


/*
 * There are four possible line/triangle intersections.
 */
typedef enum {
	tri_no_intersection,	// nothing changed
	tri_infront,		// segment is in front of the triangle
	tri_hidden,		// segment is completely occluded
	tri_clipped,		// segment is partially occluded on one end
	tri_split,		// segment is partially occluded in the middle
} tri_intersect_t;


tri_intersect_t
tri_seg_intersect(
	const tri_t * tri,
	seg_t * s,
	seg_t ** new_seg // only if tri_split
);

v3_t
tri_bary_coord(
	const tri_t * const t,
	const v3_t * const p
);

#endif
