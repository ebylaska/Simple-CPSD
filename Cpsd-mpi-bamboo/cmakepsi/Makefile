BINDIR = ../bin
SRCDIR = .
OBJECTS = makepsi.o
DOBJECTS = dmakepsi.o

LIBDIR = ../lib
LIBRARY = -lblas

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = gcc


all: single

single: ${OBJECTS}
	gcc -o cmakepsi.x ${CFLAGS} ${OBJECTS} ${LIBRARY}
	mv cmakepsi.x ${BINDIR}

double:  ${OBJECTS}
	gcc -o cmakepsi.x ${CFLAGS} ${DOBJECTS} ${LIBRARY}


clean:
	rm -f core ${OBJECTS} ${BINDIR}/cmakepsi.x
