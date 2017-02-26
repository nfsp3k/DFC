CC=gcc
CFLAGS=-O3 -g -Wall -Werror -march=native -msse4.2 --coverage

all:
	$(CC) $(CFLAGS) test.c dfc.c -o test
	gcov dfc.c test.c
clean:
	rm -rf test *~
