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
	t->min = v3_min(v3_min(t->p[0], t->p[1]), t->p[2]);
	t->max = v3_max(v3_max(t->p[0], t->p[1]), t->p[2]);

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
		if(iter->min.p[2] > t->min.p[2])
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
	fprintf(stderr, "{{{%+5.1f,%+5.1f,%5.1f }},{{%+5.1f,%+5.1f,%5.1f }},{{%+5.1f,%+5.1f,%5.1f }}}\n", // norm %.3f %.3f %.3f\n",
		t->p[0].p[0],
		t->p[0].p[1],
		t->p[0].p[2],
		t->p[1].p[0],
		t->p[1].p[1],
		t->p[1].p[2],
		t->p[2].p[0],
		t->p[2].p[1],
		t->p[2].p[2]
		//t->normal.p[0],
		//t->normal.p[1],
		//t->normal.p[2]
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
			if (!v3_eq(t0->p[i], t1->p[j]))
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
		// these are likely small triangles that can be ignored
		if (tri_debug > 3)
		{
			fprintf(stderr, "uh, three points match?\n");
			tri_print(t0);
			tri_print(t1);
		}
		return 0;
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


/*
 * Find the barycentric coordinates for point of an XY coordinate in a triangle.
 *
 * p can be written as a combination of t01 and t02,
 * p - t0 = a * (t1 - t0) + b * (t2 - t0)
 * setting t0 to 0, this becomes:
 * p = a * t1 + b * t2
 * which is two equations with two unknowns
 *
 * The x and y coordinates are based on the two sides of the triangle
 * and the z coordinate is the screen coordinate z in the triangle.
 */
v3_t
tri_bary_coord(
	const tri_t * const t,
	const v3_t * const p
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

	float a = (px * t2y - py * t2x) / (t1x * t2y - t2x * t1y);
	float b = (px * t1y - py * t1x) / (t2x * t1y - t1x * t2y);
	v3_t v = {{
		a,
		b,
		t->p[0].p[2] + a * t1z + b * t2z,
	}};

	//v.p[2] = 1.0 - v.p[0] - v.p[1];

	return v;
}

// Returns true if the bary centry point is inside the triangle
// which means that a and b are non-negative and sum to less than 1
int
tri_bary_inside(
	const v3_t p
)
{
	const float a = p.p[0];
	const float b = p.p[1];

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

	// compute the z intercept for each on the two different coordinates
	if(l0_int)
	{
		*l0_int = (v3_t){{
			p0_x + r0 * s1_x,
			p0_y + r0 * s1_y,
			p0_z + r0 * (p1_z - p0_z)
		}};
	}

	if(l1_int)
	{
		*l1_int = (v3_t){{
			p2_x + r1 * s2_x,
			p2_y + r1 * s2_y,
			p2_z + r1 * (p3_z - p2_z)
		}};
	}

	return r0;
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



/*
	tri_no_intersection,	// nothing changed
	tri_infront,		// segment is in front of the triangle
	tri_hidden,		// segment is completely occluded
	tri_clipped,		// segment is partially occluded on one end
	tri_split,		// segment is partially occluded in the middle
*/
tri_intersect_t
tri_seg_intersect(
	const tri_t * t,
	seg_t * s,
	seg_t ** new_seg // only if tri_split
)
{
	// avoid processing nearly empty segments
	const float seg_len = v3_len(&s->p[0], &s->p[1]);
	if (seg_len < EPS)
		return tri_hidden;

	const v3_t p_max = v3_max(s->p[0], s->p[1]);
	const v3_t p_min = v3_min(s->p[0], s->p[1]);

	// if the segment is closer than the triangle,
	// then we no longer have to check any further into
	// the zlist (it is sorted by depth).
	if (p_max.p[2] <= t->min.p[2])
		return tri_infront;

	// check for four quadrant outside the bounding box
	// of the triangle min/max, which would have no chance
	// of intersecting with the triangle
	if (p_min.p[0] < t->min.p[0]
	&&  p_max.p[0] < t->min.p[0])
		return tri_no_intersection;
	if (p_min.p[1] < t->min.p[1]
	&&  p_max.p[1] < t->min.p[1])
		return tri_no_intersection;
	if (p_min.p[0] > t->max.p[0]
	&&  p_max.p[0] > t->max.p[0])
		return tri_no_intersection;
	if (p_min.p[1] > t->max.p[1]
	&&  p_max.p[1] > t->max.p[1])
		return tri_no_intersection;

	// there is a possibility that this line crosses the triangle
	// compute the coordinates in triangle space
	const v3_t tp0 = tri_bary_coord(t, &s->p[0]);
	const v3_t tp1 = tri_bary_coord(t, &s->p[1]);

	// if both are inside and not both on the same edge of
	// the triangle, then the segment is totally hidden.
	if (tri_bary_inside(tp0) && tri_bary_inside(tp1))
	{
		// if the segment z is closer than the triangle z
		// then the segment is in front of the triangle
		if (s->p[0].p[2] < tp0.p[2] && s->p[1].p[2] < tp1.p[2])
			return tri_no_intersection;

		// if the barycentric coord is 0 for the same edge
		// for both points, then it is on the original line
		if (tp0.p[0] < EPS && tp1.p[0] < EPS)
			return tri_no_intersection;
		if (tp0.p[1] < EPS && tp1.p[1] < EPS)
			return tri_no_intersection;

		// compute the third barycentric coordinate and check
		float c0 = 1.0 - tp0.p[0] - tp0.p[1];
		float c1 = 1.0 - tp1.p[0] - tp1.p[1];
		if (c0 < EPS && c1 < EPS)
			return tri_no_intersection;

		// it is not on an edge and not infront of the triangle
		// so the segment is totally occluded
		return tri_hidden;
	}

	// find the intersection point for each of the three
	// sides of the triangle
	v3_t is[3] = {}; // 3d point of segment intercept
	v3_t it[3] = {}; // 3d point of triangle intercept
	float ratio[3] = {}; // length along the line
	int intersections = 0;

	for(int j = 0 ; j < 3 ; j++)
	{
		ratio[intersections] = hidden_intersect(
			&s->p[0], &s->p[1],
			&t->p[j], &t->p[(j+1)%3],
			&is[intersections],
			&it[intersections]
		);

		if (ratio[intersections] < 0)
			continue;

		// if the segment intersection is closer than the
		// triangle intersection, this does not count as
		// an intersection and we can ignore it.
		if (is[intersections].p[2] < it[intersections].p[2])
			continue;

		if (tri_debug >= 2)
		{
			fprintf(stderr, "%d ratio=%.2f %+6.1f", j, ratio[intersections], it[intersections].p[2]);
			v3_print(is[intersections]);
		}
		intersections++;
	}

	// check for duplicate intersections, which happens if
	// the lines go through at precisely the corners
	// this might mean that we hit exactly at one
	// point and two of the points are the same
	if (intersections == 3)
	{
		if (v3_eq(is[0], is[2]))
			intersections--;
		else
		if (v3_eq(is[1], is[2]))
			intersections--;
		else
		if (v3_eq(is[0], is[1]))
		{
			intersections--;
			is[1] = is[2];
			it[1] = it[2];
		}
	}

	if (intersections == 2 && v3_eq(is[0], is[1]))
		intersections--;

	// no intersections? there is nothing to do
	if (intersections == 0)
		return tri_no_intersection;

	// three intersections? maybe a very small triangle
 	if (intersections == 3)
	{
		fprintf(stderr, "THREE INTERSECTIONS?\n");
		fprintf(stderr, "is0="); v3_print(is[0]);
		fprintf(stderr, "is1="); v3_print(is[1]);
		fprintf(stderr, "is2="); v3_print(is[2]);
		svg_line("#0000FF", s->p[0].p, s->p[1].p, 8);
		return tri_no_intersection;
	}

	if (intersections == 1)
	{
		if (tri_bary_inside(tp0))
		{
			// if the intercept point on the segment is
			// closer than the intercept point on the triangle edge,
			// then there is no occlusion
			if (is[0].p[2] <= it[0].p[2])
				return tri_no_intersection;

			// clipped from intersection to p1
			s->p[0] = is[0];
			return tri_clipped;
		}

		if (tri_bary_inside(tp1))
		{
			// if the intercept point on the segment is
			// closer than the intercept point on the triangle edge,
			// then there is no occlusion
			if (is[0].p[2] < it[0].p[2])
				return tri_no_intersection;

			// clipped from p0 to intersection
			s->p[1] = is[0];
			return tri_clipped;
		}

		// something isn't right. maybe we have a small triangle?
		fprintf(stderr, "ONE INTERSECTION?");
/*
		svg_line("#FFFF00", s->p[0].p, s->p[1].p, 20);
		svg_line("#000080", t->p[0].p, t->p[1].p, 8);
		svg_line("#000080", t->p[1].p, t->p[2].p, 8);
		svg_line("#000080", t->p[2].p, t->p[0].p, 8);
*/
		seg_print(s);
		v3_print(tp0);
		v3_print(tp1);
		tri_print(t);
		*new_seg = seg_new(is[0], s->p[1]);
		s->p[1] = is[0];
		return tri_split;
	}

	// two intersections: find the one that is closer to p0
	// modify the existing segment and create a new segment
	const float d00 = v3_mag2(v3_sub(is[0], s->p[0]));
	const float d01 = v3_mag2(v3_sub(is[1], s->p[0]));
	const float d10 = v3_mag2(v3_sub(is[0], s->p[1]));
	const float d11 = v3_mag2(v3_sub(is[1], s->p[1]));

	// if any of the intersections points are zero from an
	// end point on the segment, then skip that part
	if (tri_debug > 4)
	{
		seg_print(s);
		tri_print(t);
		fprintf(stderr, "d: %f %f %f %f\n", d00, d01, d10, d11);
	}

	if (d00 < EPS && d11 < EPS)
		return tri_hidden;
	if (d01 < EPS && d10 < EPS)
		return tri_hidden;
	
	if (d00 < EPS)
	{
		s->p[0] = is[1];
		return tri_clipped;
	} else
	if (d01 < EPS)
	{
		s->p[0] = is[0];
		return tri_clipped;
	} else
	if (d10 < EPS)
	{
		s->p[1] = is[1];
		return tri_clipped;
	} else
	if (d11 < EPS)
	{
		s->p[1] = is[0];
		return tri_clipped;
	}

	// neither end points match, so create a new segment
	// that excludes the space covered by the triangle.
	// determine which is closer to point is[0]
	if (d00 < d01)
	{
		// p0 is closer to is0, so new segment is is1 to p1
		*new_seg = seg_new(is[1], s->p[1]);
		s->p[1] = is[0];
	} else {
		// p0 is closer to is1, so new segment is is0 to p1
		*new_seg = seg_new(is[0], s->p[1]);
		s->p[1] = is[1];
	}

	if (tri_debug > 3)
	{
		fprintf(stderr, "SPLIT: ");
		seg_print(*new_seg);
	}
	return tri_split;
}


int
tri_seg_hidden(
	const tri_t * zlist,
	seg_t * s,
	seg_t ** slist_visible
)
{
	int count = 0;
	if (tri_debug > 2)
	{
		fprintf(stderr, "TEST: ");
		seg_print(s);
	}

	for( const tri_t * t = zlist ; t ; t = t->next )
	{
		seg_t * new_seg = NULL;
		tri_intersect_t type = tri_seg_intersect(t, s, &new_seg);

		// if there is no intersection or if the segment has
		// been clipped on one side, keep looking
		if (type == tri_no_intersection)
			continue;
		if (type == tri_clipped)
		{
			//seg_print(s);
			continue;
		}

		// if this segment is infront of this triangle then we can
		// stop searching
		if (type == tri_infront)
			break;

		// if this segment is totally occluded, we're done
		if (type == tri_hidden)
			return count;

		// if this line has been split into two, process the
		// new segment starting at the next triangle since it
		// has already intersected this one
		if (type == tri_split)
		{
			static int recursive;
			int new_count = tri_seg_hidden(
				t->next,
				new_seg,
				slist_visible
			);
			count += new_count;
			continue;
		}

		fprintf(stderr, "unknown type %d\n", type);
		return -1;
	}

	// we've reached the end and it is still visible
	s->next = *slist_visible;
	*slist_visible = s;
	return ++count;
}
