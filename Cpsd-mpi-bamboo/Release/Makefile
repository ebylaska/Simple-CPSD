SUBDIRECTORIES = ../aimd_size ../fftpack ../fftpack3d ../sfftpack3d ../ckbppv3 ../cmakepsi ../ccpsdv3 ../ccpsdv3-mpi
EXAMPLEDIRECTORIES = ../Example0-single-serial ../Example0-double-serial ../Example1-double-serial
PARALLELDIRECTORIES = ../Example0-double-mpi

ALLDIRECTORIES = $(SUBDIRECTORIES) $(EXAMPLEDIRECTORIES) $(PARALLELDIRECTORIES)

.PHONY: all onlycode parallel $(SUBDIRECTORIES) dcpsdv3 dcpsdv3-mpi examples parallel_examples $(EXAMPLEDIRECTORIES) $(PARALLELDIRECTORIES)

all: $(SUBDIRECTORIES) dcpsdv3 dcpsdv3-mpi examples parallel_examples

parallel: $(SUBDIRECTORIES) dcpsdv3 dcpsdv3-mpi parallel_examples

onlycode: $(SUBDIRECTORIES) dcpsdv3 dcpsdv3-mpi

examples: $(EXAMPLEDIRECTORIES)

parallel_examples: $(PARALLELDIRECTORIES)

$(SUBDIRECTORIES):
	$(MAKE) -C $@ 

$(EXAMPLEDIRECTORIES):
	$(MAKE) -C $@

$(PARALLELDIRECTORIES):
	$(MAKE) -C $@

dcpsdv3:
	$(MAKE) -C ../ccpsdv3 cleanobjects
	$(MAKE) -C ../ccpsdv3 double

dcpsdv3-mpi:
	$(MAKE) -C ../ccpsdv3-mpi cleanobjects
	$(MAKE) -C ../ccpsdv3-mpi double


clean: 
	for dir in $(ALLDIRECTORIES); do \
		$(MAKE) -C $$dir clean; \
	done
