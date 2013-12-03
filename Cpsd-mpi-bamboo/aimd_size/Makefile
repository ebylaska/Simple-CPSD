BINDIR = ../bin

SOURCES = aimd_size.f

OBJECTS = aimd_size.o

FFLAGS = -c -O3

.SUFFIXES:   .f .o

.f.o:
	gfortran ${FFLAGS}  $<

all: ${OBJECTS}
	gfortran -o aimd_size ${OBJECTS} -lblas
	mv aimd_size  ${BINDIR}/aimd_size

clean:
	rm -f *.o *.a *.trace core ${BINDIR}/aimd_size



