CFLAGS += -Wall -Ofast
LDLIBS += -lm

ifdef ENABLE_OPENMP
CFLAGS += -fopenmp
endif

OBJECTS = tiff.o raytracing.o main.o

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDLIBS) $(LDFLAGS)

launch: main
	./main

format:
	clang-format -i *.c *.h

clean:
	rm *.o
