#include <stdio.h>
#include "camera.h"
#include "v3.h"

int main(void)
{
	v3_t lookat = {{0,0,0}};
	v3_t eye = {{0,-20,0}};
	v3_t up = {{0,0,1}};
	float fov = 20;

	camera_t * camera = camera_new(eye, lookat, up, fov);

	for (float x = -5 ; x <= 5 ; x += 5)
	for (float y = -25 ; y < 55 ; y += 5 )
	{
		v3_t v_in = {{ x, y, 5 }};
		v3_t v_out;
		int onscreen = camera_project(camera, &v_in, &v_out);
		if (!onscreen && y > eye.p[1])
			fprintf(stderr, "false positive\n");
		if (onscreen && y < eye.p[1])
			fprintf(stderr, "false negative\n");
	}
}
