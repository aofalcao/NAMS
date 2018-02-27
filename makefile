CC=gcc
CFLAGS = -I. -O2 -lm -lstdc++

namsmake: nams.o hungarian.o Main.o
	$(CC) -o nams nams.o hungarian.o Main.o $(CFLAGS)
