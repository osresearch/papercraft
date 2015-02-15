/** \file
 * Generate an svg file with the polygonal faces.
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
#include "stl_3d.h"

static const char * stroke_string
	= "stroke-width=\"0.1px\" fill=\"none\"";

typedef struct
{
	v3_t origin;
	v3_t x;
	v3_t y;
	v3_t z;
} refframe_t;


static void
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


static void
svg_line(
	const refframe_t * const ref,
	const v3_t p1,
	const v3_t p2
)
{
	// project p1 and p2 into the plane
	double x1, y1, x2, y2;

	v3_project(ref, p1, &x1, &y1);
	v3_project(ref, p2, &x2, &y2);
	const char * color = "#FF0000";

	printf("<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%s\" %s/>\n",
		x1, y1,
		x2, y2,
		color,
		stroke_string
	);
}


static void
svg_circle(
	const double x,
	const double y,
	const double rad,
	const char * const color
)
{
	printf("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"%s\" %s/>\n",
		x,
		y,
		rad,
		color,
		stroke_string
	);
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
inset(
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


int
main(void)
{
	stl_3d_t * const stl = stl_3d_parse(STDIN_FILENO);
	if (!stl)
		return EXIT_FAILURE;
	const double inset_distance = 8;
	const double hole_radius = 1.5;

	int * const face_used = calloc(sizeof(*face_used), stl->num_face);

	// for each vertex, find the coplanar triangles
	// \todo: do coplanar bits
	printf("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");
	printf("<g transform=\"scale(3.543307)\"><!-- scale to mm -->\n");

	const stl_vertex_t ** const vertex_list = calloc(sizeof(*vertex_list), stl->num_vertex);

	for(int i = 0 ; i < stl->num_face ; i++)
	{
		if (face_used[i])
			continue;

		const stl_face_t * const f = &stl->face[i];
		const int vertex_count = stl_trace_face(
			stl,
			f,
			vertex_list,
			face_used
		);

		fprintf(stderr, "%d: %d vertices\n", i, vertex_count);

		// generate a refernce frame based on this face
		refframe_t ref = {
			.origin = f->vertex[0]->p,
		};

		const v3_t dx = v3_norm(v3_sub(f->vertex[1]->p, ref.origin));
		const v3_t dy = v3_norm(v3_sub(f->vertex[2]->p, ref.origin));

		ref.x = dx;
		ref.z = v3_norm(v3_cross(dx, dy));
		ref.y = v3_norm(v3_cross(ref.x, ref.z));

		printf("<!-- face %d --><g>\n", i);

		// generate the polygon outline (should be one path?)
		for (int j = 0 ; j < vertex_count ; j++)
			svg_line(
				&ref,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p
			);

		// generate the inset mounting holes
		for (int j = 0 ; j < vertex_count ; j++)
		{
			double x, y;
			inset(&ref, inset_distance, &x, &y,
				vertex_list[(j+0) % vertex_count]->p,
				vertex_list[(j+1) % vertex_count]->p,
				vertex_list[(j+2) % vertex_count]->p
			);
			svg_circle(x, y, hole_radius, "#00ff00");
		}

		printf("</g>\n");
	}

	printf("</g></svg>\n");

	return 0;
}
