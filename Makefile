CFLAGS += -Wall -Ofast -MMD -MP -Iinclude
LDLIBS += -lm

ifdef ENABLE_OPENMP
CFLAGS += -fopenmp
endif

SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst src/%.c,obj/%.o,$(SOURCES))
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))

-include $(DEPENDS)  # re-compile when headers change

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDLIBS) $(LDFLAGS)

launch: main
	./main

format:
	clang-format -i src/*.c include/*.h

clean:
	rm $(OBJECTS) $(DEPENDS)
