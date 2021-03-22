CC = gcc 
LD = gcc
CFLAGS = -O3 -march=native -ffast-math -funroll-loops -g -Wall -Werror 
LDFLAGS = -lm
RM = rm -f
EXECUTABLE = galsim

galsim: galsim.c 
	$(CC) $(CFLAGS) -o $(EXECUTABLE) galsim.c $(LDFLAGS)

clean:
	$(RM) ./galsim