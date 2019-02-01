SOURCES = lulesh_mk.c
TARGET = lulesh
CFLAGS = -std=c99 -lm -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET)

clean:
	rm $(TARGET)
