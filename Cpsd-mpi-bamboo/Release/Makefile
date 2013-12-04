SUBDIRECTORIES = ../aimd_size ../fftpack ../fftpack3d ../sfftpack3d ../ckbppv3 ../cmakepsi ../ccpsdv3
EXAMPLEDIRECTORIES = ../Example0-single-serial ../Example0-double-serial

ALLDIRECTORIES = $(SUBDIRECTORIES) $(EXAMPLEDIRECTORIES)

.PHONY: all $(SUBDIRECTORIES) dcpsdv3 examples $(EXAMPLEDIRECTORIES)

all: $(SUBDIRECTORIES) dcpsdv3 examples

examples: $(EXAMPLEDIRECTORIES)

$(SUBDIRECTORIES):
	$(MAKE) -C $@ 

$(EXAMPLEDIRECTORIES):
	$(MAKE) -C $@

dcpsdv3:
	$(MAKE) -C ../ccpsdv3 cleanobjects
	$(MAKE) -C ../ccpsdv3 double

clean: 
	for dir in $(ALLDIRECTORIES); do \
		$(MAKE) -C $$dir clean; \
	done
