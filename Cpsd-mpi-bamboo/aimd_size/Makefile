BINDIR = ../bin

SOURCES = aimd_size.f
ifeq ($(HOSTTYPE),x86_64-linux)
   BLASLIB = -L/home/bylaska/lib -lblas
else
   BLASLIB = -lblas
endif
OBJECTS = aimd_size.o

FFLAGS = -c -O3

.SUFFIXES:   .f .o

.f.o:
	gfortran ${FFLAGS}  $<

all: ${OBJECTS}
	gfortran -o aimd_size ${OBJECTS} ${BLASLIB}
	mv aimd_size  ${BINDIR}/aimd_size

clean:
	rm -f *.o *.a *.trace core ${BINDIR}/aimd_size



