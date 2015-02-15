/** \file
 * STL file format.
 *
 * Parse an STL file into an easily traversed structure.
 */
#ifndef _stl3d_h_
#define _stl3d_h_

#include "v3.h"

typedef struct stl_vertex stl_vertex_t;
typedef struct stl_face stl_face_t;

#define STL_MAX_FACES 64

struct stl_vertex {
	v3_t p;
	int num_face;
	stl_face_t *face[STL_MAX_FACES];
	int face_num[STL_MAX_FACES]; // which vertex on the face
};

struct stl_face
{
	stl_vertex_t * vertex[3];
	stl_face_t * face[3];
	double angle[3];
};


typedef struct
{
	int num_vertex;
	stl_vertex_t * vertex;

	int num_face;
	stl_face_t * face;
} stl_3d_t;


stl_3d_t *
stl_3d_parse(
	int fd
);


/** Generate the list of vertices that are coplanar given a starting
 * vertex in the stl file.
 *
 * vertex_list should have enough space to contain at least as many
 * vertices as in the stl file.
 *
 * if face_used is not null it will be populated with which faces
 * have been traversed during the search.  it should have enough size
 * to contain all of the faces in the file.
 */
int
stl_trace_face(
	const stl_3d_t * const stl,
	const stl_face_t * const f_start,
	const stl_vertex_t ** vertex_list,
	int * const face_used,
	int start_vertex
);


typedef struct
{
	v3_t origin;
	v3_t x;
	v3_t y;
	v3_t z;
} refframe_t;


void
refframe_init(
	refframe_t * ref,
	const v3_t p0,
	const v3_t p1,
	const v3_t p2
);


/** Project a 3D point onto a 2D space */
void
v3_project(
	const refframe_t * const ref,
	const v3_t p_in,
	double * const x_out,
	double * const y_out
);


#endif
