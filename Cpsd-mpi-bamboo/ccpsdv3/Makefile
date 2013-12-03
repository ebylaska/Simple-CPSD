BINDIR = ../bin
SRCDIR = .
OBJECTS = main.o cpsd.o dirac_exchange.o vosko.o message.o get_cube.o get_g.o get_masker.o strfac.o gcdotc.o gsdot.o gssum.o eigsrt.o second.o

LIBDIR = ../lib
LIBRARYd = -lfftpack3d -lfftpack -llapack -lblas
LIBRARYs = -lfftpack3ds -lfftpack -llapack -lblas

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = gcc


all: single

single: setsingle ${OBJECTS}
	gcc -o ccpsdv3.x ${CFLAGS} ${OBJECTS} ${LIBRARYs}
	mv ccpsdv3.x ${BINDIR}

double: setdouble ${OBJECTS}
	gcc -o ccpsdv3.x ${CFLAGS} ${OBJECTS} ${LIBRARYd}

cpsd.x: ${OBJECTS}
	gcc -o ccpsdv3.x ${CFLAGS} ${OBJECTS} ${LIBRARY}

setsingle:
	cp float_single.h float.h

setdouble:
	cp float_double.h float.h

clean:
	rm -f core ${OBJECTS} ${BINDIR}/ccpsdv3.x
