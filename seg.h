#ifndef _seg_h_
#define _seg_h_

#include "v3.h"

typedef struct _seg_t seg_t;
struct _seg_t {
	v3_t p[2];
	v3_t src[2];
	seg_t * next;
};

static inline seg_t *
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


static inline void
seg_print(
	const seg_t * const s
)
{
	fprintf(stderr, "%.0f,%.0f -> %.0f,%.0f (was %.0f,%.0f -> %.0f,%.0f)\n",
		s->p[0].p[0],
		s->p[0].p[1],
		s->p[1].p[0],
		s->p[1].p[1],
		s->src[0].p[0],
		s->src[0].p[1],
		s->src[1].p[0],
		s->src[1].p[1]
	);
}


#endif
