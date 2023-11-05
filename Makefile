CFLAGS += -Wall -Ofast -MMD -MP -I./include
LDLIBS += -lm

ifdef ENABLE_OPENMP
CFLAGS += -fopenmp
endif

SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst %.c,%.o,$(SOURCES))
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))

-include $(DEPENDS)  # re-compile when headers change

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDLIBS) $(LDFLAGS)

launch: main
	./main

format:
	clang-format -i src/*.c include/*.h

clean:
	rm $(OBJECTS) $(DEPENDS)
