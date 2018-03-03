#ifndef _svg_h_
#define _svg_h_

static inline void
svg_line(
	const char * color,
	const float * p1,
	const float * p2,
	float thick
)
{
	// invert the sense of y
	printf("<line x1=\"%fpx\" y1=\"%fpx\" x2=\"%fpx\" y2=\"%fpx\" stroke=\"%s\" stroke-width=\"%.1fpx\"/>\n",
		p1[0],
		-p1[1],
		p2[0],
		-p2[1],
		color,
		thick
	);
}

static inline void
svg_circle(
	const char * color,
	const float  p1,
	const float  p2,
	float radius
)
{
	// invert the sense of y
	printf("<circle x=\"%fpx\" y=\"%fpx\" radius=\"%fpx\" stroke=\"%s\" stroke-width=\"%.1fpx\"/>\n",
		p1,
		-p2,
		radius,
		color,
		1.0
	);
}

#endif
