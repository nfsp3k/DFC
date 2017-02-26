CC=gcc
CFLAGS=-O3 -g -Wall -Werror -march=native -msse4.2

all:
	$(CC) $(CFLAGS) test.c dfc.c -o test

clean:
	rm -rf test *~
