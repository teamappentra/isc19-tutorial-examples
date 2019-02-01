SOURCES = atmux.c lib/Matrix2D.c lib/Vector.c lib/CRSMatrix.c
TARGET = atmux
CFLAGS = -std=c99 -O2 -Ilib

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET) 10000

clean:
	rm $(TARGET)
