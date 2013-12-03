LIBDIR = ../lib

SOURCES = fftpack3d.f

OBJECTS = fftpack3d.o

FFLAGS = -c -O3

.SUFFIXES:   .f .o

.f.o:
	gfortran ${FFLAGS}  $<

all: ${OBJECTS}
	ar rc libfftpack3d.a  ${OBJECTS}
	mv libfftpack3d.a  ${LIBDIR}

clean:
	rm -f *.o *.a *.trace core ${LIBDIR}/libfftpack3d.a



