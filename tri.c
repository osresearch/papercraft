/*
 * Triangle manipulations
 */
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "tri.h"


int tri_debug = 0;


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
	t->prev = zlist;
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


// Compute the 2D area of a triangle in screen space
// using Heron's formula
float
tri_area_2d(
	const tri_t * const t
)
{
	const float a = v3_dist_2d(&t->p[0], &t->p[1]);
	const float b = v3_dist_2d(&t->p[1], &t->p[2]);
	const float c = v3_dist_2d(&t->p[2], &t->p[0]);
	const float s = (a + b + c) / 2;

	return sqrt(s * (s-a) * (s-b) * (s-c));
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
 *
 * Returns true if the point is inside the triangle
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

//fprintf(stderr, "%d: processing segment ", recursive); seg_print(s);
	fprintf(stderr, "--- recursive %d\n", recursive);
	seg_print(s);

	for( const tri_t * t = zlist ; t ; t = t->next )
	{
		// if the segment is closer than the triangle,
		// then we no longer have to check any further into
		// the zlist (it is sorted by depth).
		if (seg_max_z <= t->min[2])
			break;

#if 0
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
#endif

		if (tri_debug >= 2)
			tri_print(t);
/*
		// if the segment is co-linear to any of the
		// triangle edges, include it
		for(int i = 0 ; i < 3 ; i++)
		{
			if (parallel(
				&s->p[0], &s->p[1],
				&t->p[i], &t->p[(i+1)%3]
			))
				goto next_segment;
		}
*/

		float z0, z1;
		int inside0 = tri_find_z(t, &s->p[0], &z0);
		int inside1 = tri_find_z(t, &s->p[1], &z1);

		if (tri_debug >= 2 && (inside0 || inside1))
		{
			fprintf(stderr, "inside %d %d\n", inside0, inside1);
		}

		// if both are inside but the segment is infront of the
		// triangle, then we retain the segment.
		// otherwies we discard the segment
		if (inside0 && inside1)
		{
			if (s->p[0].p[2] <= z0
			&&  s->p[1].p[2] <= z1)
				continue;
			if (tri_debug >= 2)
				fprintf(stderr, "BOTH INSIDE\n");
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

			if (tri_debug >= 2)
				fprintf(stderr, "%d ratio=%.2f\n", j, ratio);
			intersections++;
		}


		// if none of them intersect, we keep looking
		if (intersections == 0)
			continue;

		if (tri_debug >= 2)
			fprintf(stderr, "%d intersections\n", intersections);

		if (intersections == 3)
		{
			// this likely means that the triangle is very, very
			// small, so let's just ignore this triangle
			if (tri_debug >= 2)
				fprintf(stderr, "Three intersections\n");
			continue;
		}


		if (intersections == 2)
		{
			// figure out how far it is to each of the intersections
			const float d00 = v3_len(&s->p[0], &is[0]);
			const float d01 = v3_len(&s->p[0], &is[1]);
			const float d10 = v3_len(&s->p[1], &is[0]);
			const float d11 = v3_len(&s->p[1], &is[1]);

			if (tri_debug >= 2)
				fprintf(stderr, "Two intersections\n");

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

next_segment:
		continue;
	}

	// if we've reached here the segment is visible
	// and should be added to the visible list
	s->next = *slist_visible;
	*slist_visible = s;
	recursive--;
}


/*
 * Fast check to see if t2 is entire occluded by t.
 */
int
tri_behind(
	const tri_t * const t,
	const tri_t * const t2
)
{
	float z0, z1, z2;
	int inside0 = tri_find_z(t, &t2->p[0], &z0);
	int inside1 = tri_find_z(t, &t2->p[1], &z1);
	int inside2 = tri_find_z(t, &t2->p[2], &z2);

	// easy check -- if none of the points are inside,
	// t2 is not entirely occluded
	if (!inside0 || !inside1 || !inside2)
		return 0;

	// are all of the intersection points ahead of t2?
	int behind0 = t2->p[0].p[2] >= z0;
	int behind1 = t2->p[1].p[2] >= z1;
	int behind2 = t2->p[2].p[2] >= z2;
	if (behind0 && behind1 && behind2)
		return 1;

	// it is a STL violation if they are not all on the
	// same side (this would indicate that t and t2 intersect
	// go ahead and prune since it will cause problems
	if (behind0 || behind1 || behind2)
	{
/*
		fprintf(stderr, "WARNING: triangles intersect %.0f %.0f %.0f inside %d %d %d behind %d %d %d\n", z0, z1, z2, inside0, inside1, inside2, behind0, behind1, behind2);
		tri_print(t);
		tri_print(t2);
*/
		return 1;
	}

	// they are all on the same side
	return 0;
}
