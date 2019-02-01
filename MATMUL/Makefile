SOURCES = matmul.c
TARGET = matmul
CFLAGS = -std=c99 -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 1000

clean:
	rm $(TARGET)
