CC = gcc
CFLAGS = -O3 -march=native -lrt
MPM_PATH=mpm_engines
DFC_FILES=../src/dfc.c ../src/dfc.h ../src/dfc_framework.h
AC_FILES=${MPM_PATH}/ac/acsmx2.c ${MPM_PATH}/ac/acsmx2.h ${MPM_PATH}/ac/util_str.c ${MPM_PATH}/ac/util_str.h

all: test
test: test.c entry.h ${DFC_FILES} ${AC_FILES} match_func.h match_func.c
	${CC} ${CFLAGS} -g -o test test.c ${DFC_FILES} ${AC_FILES} match_func.c

clean:
	rm -rf test
