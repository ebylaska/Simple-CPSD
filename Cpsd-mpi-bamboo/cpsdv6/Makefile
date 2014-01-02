BINDIR = ../bin

OBJECTS = main.o cpsd.o dirac_exchange.o vosko.o message.o get_cube.o eigsrt.o second.o

ifeq ($(HOSTTYPE),x86_64-linux)
   BLASLIB = -L/home/bylaska/lib -llapack -lblas -lm -lgfortran
else
   BLASLIB = -llapack -lblas
endif


LIBDIR = ../lib
LIBRARYs = -lfftpack3ds -lfftpack ${BLASLIB} 
LIBRARYd = -lfftpack3d -lfftpack ${BLASLIB}

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = gcc


all: single

single: setsingle ${OBJECTS}
	gcc -o scpsdv6.x ${CFLAGS} ${OBJECTS} ${LIBRARYs}
	mv scpsdv6.x ${BINDIR}

double: setdouble ${OBJECTS}
	gcc -o dcpsdv6.x ${CFLAGS} ${OBJECTS} ${LIBRARYd}
	mv dcpsdv6.x ${BINDIR}


setsingle:
	cp float_single.h float.h

setdouble:
	cp float_double.h float.h

cleanobjects:
	rm -f ${OBJECTS} 

clean:
	rm -f core ${OBJECTS} ${BINDIR}/scpsdv6.x ${BINDIR}/dcpsdv6.x
