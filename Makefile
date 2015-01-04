CFLAGS = \
	-std=c99 \
	-g \
	-O3 \
	-W \
	-Wall \

LDLIBS = \
	-lm \

all: unfold

unfold: unfold.o

clean:
	$(RM) *.o
