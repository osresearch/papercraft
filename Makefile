CFLAGS = \
	-std=c99 \
	-g \
	-O3 \
	-W \
	-Wall \
	-Wp,-MMD,$(dir $@).$(notdir $@).d \
	-Wp,-MT,$@ \

LDLIBS = \
	-lm \

all: unfold wireframe corners

unfold: unfold.o
wireframe: wireframe.o
corners: corners.o stl_3d.o

clean:
	$(RM) *.o

-include .*.o.d
