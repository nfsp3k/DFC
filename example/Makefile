CC=gcc
CFLAGS= -g -O3 -msse4 -Wall -march=native -lrt
DFC_PATH=../src


all: test

test: test.c Makefile ${DFC_PATH}/dfc.c ${DFC_PATH}/dfc.h ${DFC_PATH}/dfc_framework.h
		${CC} ${CFLAGS} -o test ${DFC_PATH}/dfc.c test.c

clean:
		rm test
