BINDIR = ../bin


OBJECTS = kbppv3.o
DOBJECTS = dkbppv3.o
ifeq ($(HOSTTYPE),x86_64-linux)
   BLASLIB = -L/home/bylaska/lib -lblas -lm -gfortran
else
   BLASLIB = -lblas
endif


LIBDIR = ../lib
LIBRARYd = $(BLASLIB)
LIBRARYs = $(BLASLIB)

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = gcc


all: single double

single: ${OBJECTS}
	gcc -o skbppv3.x ${CFLAGS} ${OBJECTS} ${LIBRARYs}
	mv skbppv3.x ${BINDIR}

double:  ${DOBJECTS}
	gcc -o dkbppv3.x ${CFLAGS} ${DOBJECTS} ${LIBRARYd}
	mv dkbppv3.x ${BINDIR}


clean:
	rm -f core ${OBJECTS} ${DOBJECTS} ${BINDIR}/skbppv3.x ${BINDIR}/dkbppv3.x
