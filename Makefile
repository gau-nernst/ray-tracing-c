CFLAGS += -Wall -Ofast -MMD -MP
LDLIBS += -lm

ifdef ENABLE_OPENMP
CFLAGS += -fopenmp
endif

OBJECTS = tiff.o pcg32.o vec3.o raytracing.o main.o
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))

-include $(DEPENDS)  # re-compile when headers change

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDLIBS) $(LDFLAGS)

launch: main
	./main

format:
	clang-format -i *.c *.h

clean:
	rm *.o *.d
