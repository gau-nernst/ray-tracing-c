CFLAGS += -std=c11 -Wall -Ofast -MMD -MP -Iinclude -Istb
LDLIBS += -lm

ifdef ENABLE_OPENMP
CFLAGS += -fopenmp
endif

SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst src/%.c,obj/%.o,$(SOURCES))
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))

-include $(DEPENDS)  # re-compile when headers change

$(shell mkdir -p obj)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDLIBS) $(LDFLAGS)

launch%: main
	./main $(patsubst launch%,%,$@)

format:
	clang-format -i src/*.c include/*.h

clean:
	rm $(OBJECTS) $(DEPENDS)
