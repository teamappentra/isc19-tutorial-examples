SOURCES = laplace.c
TARGET = laplace
CFLAGS = -std=c99 -lm -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 200

clean:
	rm $(TARGET)
