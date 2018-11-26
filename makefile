CC=gcc
CFLAGS = -I. -O2 -lm -lstdc++

namsmake: nams.o hungarian.o Main.o awts.o
	$(CC) -o nams nams.o awts.o hungarian.o Main.o $(CFLAGS)
