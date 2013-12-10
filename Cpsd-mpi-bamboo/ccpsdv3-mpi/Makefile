BINDIR = ../bin

OBJECTS = main.o cpsd.o dirac_exchange.o vosko.o message.o get_cube.o get_g.o get_masker.o strfac.o gcdotc.o gsdot.o gssum.o eigsrt.o second.o d3db-mpi.o hilbert.o hcurve.o olist.o Parallel.o

LIBDIR = ../lib
LIBRARYs = -lfftpack3ds -lfftpack -llapack -lblas
LIBRARYd = -lfftpack3d -lfftpack -llapack -lblas

OPTIONS = -O3  -L$(LIBDIR)
CFLAGS =  $(OPTIONS) 
CC = mpicc


all: single

single: setsingle ${OBJECTS}
	mpicc -o scpsdv3-mpi.x ${CFLAGS} ${OBJECTS} ${LIBRARYs}
	mv scpsdv3-mpi.x ${BINDIR}

double: setdouble ${OBJECTS}
	mpicc -o dcpsdv3-mpi.x ${CFLAGS} ${OBJECTS} ${LIBRARYd}
	mv dcpsdv3-mpi.x ${BINDIR}


setsingle:
	cp float_single.h float.h

setdouble:
	cp float_double.h float.h

cleanobjects:
	rm -f ${OBJECTS} 

clean:
	rm -f core ${OBJECTS} ${BINDIR}/scpsdv3.x ${BINDIR}/dcpsdv3.x
