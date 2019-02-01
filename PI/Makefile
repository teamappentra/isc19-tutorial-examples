SOURCES = pi.c
TARGET = pi
CFLAGS = -std=c99 -lm -O2

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 100000000

clean:
	rm $(TARGET)
