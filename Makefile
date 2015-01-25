CFLAGS = \
	-std=c99 \
	-g \
	-O3 \
	-W \
	-Wall \

LDLIBS = \
	-lm \

all: unfold wireframe

unfold: unfold.o
wireframe: wireframe.o

clean:
	$(RM) *.o
