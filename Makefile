CFLAGS += -Wall -O3 -Ofast
LDLIBS += -lm

OBJECTS = tiff.o raytracing.o main.o

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDLIBS)

launch: main
	./main

format:
	clang-format -i *.c *.h
