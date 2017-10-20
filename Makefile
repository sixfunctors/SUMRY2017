# Sets the default rules for make
CC=g++
CFLAGS=-g3 -Wall -std=c99 -pedantic


Gold: bigint.o usefcns.o
	${CC} ${CFLAGS} -o Gold ThresholdEnumeration.cpp bigint.o usefcns.o

bigint.o:
	${CC} ${CFLAGS} -c bigint.cpp

usefcns.o:
	${CC} ${CFLAGS} -c usefcns.cpp

clean:
	rm -f *.o