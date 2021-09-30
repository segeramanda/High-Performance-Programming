CC = gcc-11 
LD = gcc
CFLAGS = -O3 -march=native -ffast-math 
LDFLAGS = -funroll-loops -lm -fopenmp
RM = rm -f
EXECUTABLE = golmp

g: golnewmp.c 
	$(CC) $(CFLAGS) -o $(EXECUTABLE) golnewmp.c $(LDFLAGS)

clean:
	$(RM) ./golmp