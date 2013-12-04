BINDIR = ../bin

OBJECTS = makepsi.o
DOBJECTS = dmakepsi.o

LIBDIR = ../lib
LIBRARY = -lblas

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = gcc


all: single double

single: ${OBJECTS}
	gcc -o smakepsi.x ${CFLAGS} ${OBJECTS} ${LIBRARY}
	mv smakepsi.x ${BINDIR}

double:  ${DOBJECTS}
	gcc -o dmakepsi.x ${CFLAGS} ${DOBJECTS} ${LIBRARY}
	mv dmakepsi.x ${BINDIR}


clean:
	rm -f core ${OBJECTS} ${DOBJECTS} ${BINDIR}/smakepsi.x ${BINDIR}/dmakepsi.x
