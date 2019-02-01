SOURCES = daxpy.c
TARGET = daxpy
CFLAGS = -std=c99 -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 100000000

clean:
	rm $(TARGET)
