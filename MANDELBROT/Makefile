SOURCES = mandelbrot.c
TARGET = mandelbrot
CFLAGS = -std=c99 -lm -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 5000

clean:
	rm $(TARGET)
