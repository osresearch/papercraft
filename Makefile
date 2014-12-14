CFLAGS = \
	-g \
	-O3 \
	-W \
	-Wall \

LDFLAGS = \
	-lm \

all: unfold

unfold: unfold.o

clean:
	$(RM) *.o
