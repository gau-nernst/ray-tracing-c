CFLAGS += -Wall

OBJECTS = tiff.o main.o

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

launch: main
	./main

format:
	clang-format -i *.c *.h
