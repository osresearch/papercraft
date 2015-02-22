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
	
	//int check = -EPS < dot && dot < +EPS;
	int check = -1 < dot && dot < +1;

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


/** Starting at a point, trace the coplanar polygon and return a
 * list of vertices.
 */
int
stl_trace_face(
	const stl_3d_t * const stl,
	const stl_face_t * const f_start,
	const stl_vertex_t ** vertex_list,
	int * const face_used,
	const int start_vertex
)
{
	const stl_face_t * f = f_start;
	int i = start_vertex;
	int vertex_count = 0;

	do {
		const stl_vertex_t * const v1 = f->vertex[(i+0) % 3];
		const stl_vertex_t * const v2 = f->vertex[(i+1) % 3];
		const stl_face_t * const f_next = f->face[i];

		fprintf(stderr, "%p %d: %f,%f,%f\n", f, i, v1->p.p[0], v1->p.p[1], v1->p.p[2]);
		if (face_used)
			face_used[f - stl->face] = 1;

		if (!f_next || f->angle[i] != 0)
		{
			// not coplanar or no connection.
			// add the NEXT vertex on this face and continue
			vertex_list[vertex_count++] = v2;
			i = (i+1) % 3;
			continue;
		}

		// coplanar; figure out which vertex on the next
		// face we start at
		int i_next = -1;
		for (int j = 0 ; j < 3 ; j++)
		{
			if (f_next->vertex[j] != v1)
				continue;
			i_next = j;
			break;
		}

		if (i_next == -1)
			abort();
		
		// move to the new face
		f = f_next;
		i = i_next;

		// keep going until we reach our starting face
		// at vertex 0.
	} while (f != f_start || i != start_vertex);

	return vertex_count;
}


void
refframe_init(
	refframe_t * ref,
	const v3_t p0,
	const v3_t p1,
	const v3_t p2
)
{
	ref->origin = p0;

	const v3_t dx = v3_norm(v3_sub(p1, ref->origin));
	const v3_t dy = v3_norm(v3_sub(p2, ref->origin));

	ref->x = dx;
	ref->z = v3_norm(v3_cross(dx, dy));
	ref->y = v3_norm(v3_cross(ref->x, ref->z));
}


void
v3_project(
	const refframe_t * const ref,
	const v3_t p_in,
	double * const x_out,
	double * const y_out
)
{
	v3_t p = v3_sub(p_in, ref->origin);

	double x = ref->x.p[0]*p.p[0] + ref->x.p[1]*p.p[1] + ref->x.p[2]*p.p[2];
	double y = ref->y.p[0]*p.p[0] + ref->y.p[1]*p.p[1] + ref->y.p[2]*p.p[2];
	double z = ref->z.p[0]*p.p[0] + ref->z.p[1]*p.p[1] + ref->z.p[2]*p.p[2];

	fprintf(stderr, "%f,%f,%f\n", x, y, z);

	*x_out = x;
	*y_out = y;
}


//  Determines the intersection point of the line defined by points A and B with the
//  line defined by points C and D.
//
//  Returns YES if the intersection point was found, and stores that point in X,Y.
//  Returns NO if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

static int
line_intersect(
	double Ax, double Ay,
	double Bx, double By,
	double Cx, double Cy,
	double Dx, double Dy,
	double *X, double *Y
)
{
	//  Fail if either line is undefined.
	if ((Ax==Bx && Ay==By) || (Cx==Dx && Cy==Dy))
		return 0;

	//  (1) Translate the system so that point A is on the origin.
	Bx-=Ax; By-=Ay;
	Cx-=Ax; Cy-=Ay;
	Dx-=Ax; Dy-=Ay;

	//  Discover the length of segment A-B.
	const double distAB=sqrt(Bx*Bx+By*By);

	//  (2) Rotate the system so that point B is on the positive X axis.
	const double theCos=Bx/distAB;
  	const double theSin=By/distAB;

  	double newX=Cx*theCos+Cy*theSin;
  	Cy  =Cy*theCos-Cx*theSin; Cx=newX;
  	newX=Dx*theCos+Dy*theSin;
  	Dy  =Dy*theCos-Dx*theSin; Dx=newX;

	//  Fail if the lines are parallel.
	if (Cy==Dy) return 0;

	//  (3) Discover the position of the intersection point along line A-B.
	const double ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

	//  (4) Apply the discovered position to line A-B in the original coordinate system.
	*X=Ax+ABpos*theCos;
	*Y=Ay+ABpos*theSin;

	return 1;
}



/** Compute the inset coordinate.
 * http://alienryderflex.com/polygon_inset/
//  Given the sequentially connected points (a,b), (c,d), and (e,f), this
//  function returns, in (C,D), a bevel-inset replacement for point (c,d).
//
//  Note:  If vectors (a,b)->(c,d) and (c,d)->(e,f) are exactly 180° opposed,
//         or if either segment is zero-length, this function will do
//         nothing; i.e. point (C,D) will not be set.

 */
void
refframe_inset(
	const refframe_t * const ref,
	const double inset_dist,
	double * const x_out,
	double * const y_out,
	const v3_t p0, // previous point
	const v3_t p1, // current point to inset
	const v3_t p2  // next point
)
{
	double a, b, c, d, e, f;
	v3_project(ref, p0, &a, &b);
	v3_project(ref, p1, &c, &d);
	v3_project(ref, p2, &e, &f);

	double c1 = c;
	double d1 = d;
	double c2 = c;
	double d2 = d;

	//  Calculate length of line segments.
	const double dx1 = c-a;
	const double dy1 = d-b;
	const double dist1 = sqrt(dx1*dx1+dy1*dy1);
	const double dx2 = e-c;
	const double dy2 = f-d;
	const double dist2 = sqrt(dx2*dx2+dy2*dy2);

	//  Exit if either segment is zero-length.
	if (dist1==0. || dist2==0.)
	{
		*x_out = *y_out = 0;
		fprintf(stderr, "inset fail\n");
		return;
	}

	//  Inset each of the two line segments.
	double insetX, insetY;

	insetX =  dy1/dist1 * inset_dist;
	a+=insetX;
	c1+=insetX;

	insetY = -dx1/dist1 * inset_dist;
	b+=insetY;
	d1+=insetY;

	insetX = dy2/dist2 * inset_dist;
	e+=insetX;
	c2+=insetX;

	insetY = -dx2/dist2 * inset_dist;
	f+=insetY;
	d2+=insetY;

	//  If inset segments connect perfectly, return the connection point.
	if (c1==c2 && d1==d2)
	{
		*x_out = c1;
		*y_out = d1;
		return;
	}

	//  Return the intersection point of the two inset segments (if any).
	if (line_intersect(a,b,c1,d1,c2,d2,e,f, x_out, y_out))
		return;

	*x_out = *y_out = 0;
	fprintf(stderr, "inset failed 2\n");
}
