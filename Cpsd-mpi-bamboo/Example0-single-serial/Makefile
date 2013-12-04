
all: C01.vpp c2.elc00 c2.msg01

C01.vpp:
	echo 2 | ../bin/skbppv3 C01

c2.elc00:
	../bin/smakepsi c2.elc00

c2.msg01:
	./cpsdv3.job

clean:
	rm -f IONIN IONOUT ELCIN ELCOUT MESSAGE CONTROL PSI0.dat C01.vpp c2.elc00 c2.elc01 c2.ion01 c2.msg01



