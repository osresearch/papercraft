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

#endif
