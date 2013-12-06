BINDIR = ../bin

OBJECTS = makepsi.o
DOBJECTS = dmakepsi.o

ifeq ($(HOSTTYPE),x86_64-linux)
   BLASLIB = -L/home/bylaska/lib -lblas -lm -lgfortran
else
   BLASLIB = -lblas
endif


LIBDIR = ../lib
LIBRARY = ${BLASLIB}

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
