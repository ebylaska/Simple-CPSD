LIBDIR = ../lib

SOURCES = fftpack3d.f
OBJECTS = fftpack3d.o

FFLAGS = -c -O3

.SUFFIXES:   .f .o

.f.o:
	gfortran ${FFLAGS}  $<

all:	${OBJECTS}
			ar rc libfftpack3ds.a  ${OBJECTS}
			mv libfftpack3ds.a  ${LIBDIR}

clean:
	rm -f *.o *.a *.trace core ${LIBDIR}/libfftpack3ds.a



