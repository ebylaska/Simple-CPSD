
all: C01.vpp c20.elc00 c20.msg01

C01.vpp:
	echo 2 | ../bin/dkbppv3 C01

c20.elc00:
	../bin/dmakepsi c20.elc00

c20.msg01:
	./cpsdv3.job

clean:
	rm -f IONIN IONOUT ELCIN ELCOUT MESSAGE CONTROL PSI0.dat C01.vpp c20.elc00 c20.elc01 c20.ion01 c20.msg01



