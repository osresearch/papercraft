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

all: unfold wireframe corners faces hiddenwire

unfold: unfold.o
wireframe: wireframe.o
corners: corners.o stl_3d.o
faces: faces.o stl_3d.o
hiddenwire: hiddenwire.o camera.o tri.o

# Test the triangle intersection code
test-intersect: test-intersect.o tri.o

clean:
	$(RM) *.o

-include .*.o.d
