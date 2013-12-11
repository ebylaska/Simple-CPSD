
all: C01.vpp c2.elc00 c2.msg01 c2.msg02 c2.msg04 c2.msg08

C01.vpp:
	echo 2 | ../bin/dkbppv3 C01

c2.elc00:
	../bin/dmakepsi c2.elc00

c2.msg01:
	./cpsdv3.job

c2.msg02:
	./cpsdv3.job2

c2.msg04:
	./cpsdv3.job4

c2.msg08:
	./cpsdv3.job8

clean:
	rm -f IONIN IONOUT ELCIN ELCOUT MESSAGE CONTROL PSI0.dat C01.vpp c2.elc00 c2.elc01 c2.ion01 c2.msg01 c2.msg02 c2.msg04 c2.msg08



