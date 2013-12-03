      PROGRAM MAIN
      IMPLICIT double precision (a-h,o-z)
      implicit integer (i-n)

      INCLUDE 'aimd.h'

*     secondary parameters that determine array size
*     NEMAX  : maximum number of electrons
*     NFFT3D : size of array for fast fourier transform
*     NSHL3D : number of lattice to be summed in the Ewald summation
      PARAMETER (NEMAX=NELC1+NELC2)
      PARAMETER (NFFT3D=(NFFT/2+1)*NFFT*NFFT)
      PARAMETER (NSHL3D=(2*NSHELL+1)**3)

      LOGICAL MOVE
      double precision LMD
      complex*16 C2,C1
      CHARACTER*2 ATOM(NKATMX)
      CHARACTER*8 FATOM(NKATMX)

      DIMENSION KATM(NATMX),VL(NWAVE,NKATMX),VNL(NWAVE,9,NKATMX)
      DIMENSION ZV(NKATMX),LMMAX(NKATMX),AMASS(NKATMX),VNLNRM(9,NKATMX)
      DIMENSION NATOM(NKATMX),LMAX(NKATMX),LMAX0(NKATMX)
      DIMENSION INDX(NWAVE)
      DIMENSION VC(NWAVE),VG(NWAVE),TG(NWAVE),G(NWAVE,3)
      DIMENSION UNITA(3,3),UNITG(3,3),RCELL(NSHL3D,3)
      DIMENSION NE(2),N1(2),N2(2),E(10),EN(2)
      DIMENSION R2(NATMX,3),R1(NATMX,3)
      DIMENSION DN(2*NFFT3D,2),XCP(2*NFFT3D,2),XCE(2*NFFT3D,2)
      DIMENSION LMD(NELC1,NELC1,2),HML(NELC1,NELC1,2),EIG(NELC1,2)
      DIMENSION XTMP(2*NFFT3D)
      DIMENSION C2(NWAVE,NEMAX),C1(NWAVE,NEMAX)
      dimension rc(0:4,nkatmx)

      COMMON            c2,c1,lmd,hml,fmass,ne
      COMMON / ION    / r2,r1,amass,ni
      COMMON / DENSTY / dn,xcp,xce,ispin
      COMMON / FOURIE / tg,vc,vg,g
      COMMON / LATTIC / unita,unitg,unit,omega,icube
      COMMON / EWALD  / rcell,cewald,rcut,ncut
      COMMON / PSEUDO / vl,vnl,vnlnrm,zv,nkatm,katm,lmmax
      COMMON / INDEX2  / indx,nida,nidb
      COMMON / WORKSP / xtmp

      DATA NATOM /NKATMX*0/

*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      CALL current_second(CPU1)
      OPEN(UNIT=10,FILE='MESSAGE',STATUS='NEW',FORM='FORMATTED')

      WRITE(10,1000)
      WRITE(10,1010)
      WRITE(10,1020)
      WRITE(10,1010)
      WRITE(10,1030)
      WRITE(10,1010)
      WRITE(10,1035)
      WRITE(10,1010)
      WRITE(10,1040)
      WRITE(10,1010)
      WRITE(10,1000)
      CALL MESSAGE(1)
      IF(NELC1.LT.NELC2) THEN
        IERR=100
        GO TO 9000
      ENDIF

*::::::::::::::::::  read parameters and flags  :::::::::::::::::::::::
      OPEN(UNIT=11,FILE='CONTROL',STATUS='OLD',FORM='FORMATTED')
      READ(11,*,ERR=9110,END=9111) ISPIN,IMOVE
      READ(11,*,ERR=9110,END=9111) ICUBE,UNIT
      READ(11,*,ERR=9110,END=9111) DT,FMASS
      READ(11,*,ERR=9110,END=9111) ITEST,ISTOP
      READ(11,*,ERR=9110,END=9111) TOLE,TOLC,TOLR
      IF(ISPIN.LT.1.OR.ISPIN.GT.2) THEN
        IERR=101
        GO TO 9000
      ENDIF
      IF(ICUBE.NE.JCUBE) THEN
        IERR=102
        GO TO 9000
      ENDIF
      MOVE=IMOVE.EQ.1
      RCUT=0.0d0
      NCUT=0


*::::::::::::::::::::::  read ionic structure  ::::::::::::::::::::::::
      OPEN(UNIT=12,FILE='IONIN',STATUS='OLD',FORM='FORMATTED')
      READ(12,*,ERR=9120,END=9121) NKATM
      IF(NKATM.GT.NKATMX) THEN
        IERR=103
        GO TO 9000
      ENDIF
      READ(12,*,ERR=9120,END=9121) (FATOM(IA),LMAX(IA),IA=1,NKATM)
      READ(12,*,ERR=9120,END=9121) NI
      IF(NI.GT.NATMX) THEN
        IERR=104
        GO TO 9000
      ENDIF

*     In order to maintain compatibility of file structures with      *
*     ab initio molecular dynamics structure,  velocity is also       *
*     read.                                                           *
      READ(12,*,ERR=9120,END=9121) (KATM(J),(R2(J,I),I=1,3),
     &                                      (R1(J,I),I=1,3),J=1,NI)
      DO 100 I=1,NI
        IF(KATM(I).GT.NKATM) THEN
          IERR=105
          GO TO 9000
        ENDIF
        NATOM(KATM(I))=NATOM(KATM(I))+1
  100 CONTINUE

*::::::::::::::::::   read pseudopotential  :::::::::::::::::::::::::::
      DO 120 IA=1,NKATM
        OPEN(UNIT=13,FILE=FATOM(IA),STATUS='OLD',FORM='UNFORMATTED')
        READ(13,ERR=9130,END=9131) ICUBE,NFFT0,NWAVE0,UNIT0
        IF(ICUBE.NE.JCUBE) THEN
          IERR=106
          GO TO 9000
        ENDIF
        IF(NFFT0.NE.NFFT) THEN
          IERR=107
          GO TO 9000
        ENDIF
        IF(NWAVE0.NE.NWAVE) THEN
          IERR=108
          GO TO 9000
        ENDIF
        IF(UNIT0.NE.UNIT) THEN
          IERR=109
          GO TO 9000
        ENDIF
*       READ(13,ERR=9130,END=9131) ATOM(IA),AMASS(IA),ZV(IA),LMAX0(IA)
        read(13,err=9130,end=9131) atom(ia),amass(ia),zv(ia),lmax0(ia)
        read(13,err=9130,end=9131) (rc(i,ia),i=0,lmax0(ia))
        IF(LMAX0(IA).NE.LMAX(IA)) THEN
          IERR=116
          GO TO 9000
        ENDIF
        LMMAX(IA)=LMAX(IA)**2
        IF(LMMAX(IA).GT.9) THEN
          IERR=117
          GO TO 9000
        ENDIF
        READ(13,ERR=9130,END=9131) (VNLNRM(L,IA),L=1,LMMAX(IA))
        READ(13,ERR=9130,END=9131) (VL(K,IA),K=1,NWAVE)
        DO 110 L=1,LMMAX(IA)
          READ(13,ERR=9130,END=9131) (VNL(K,L,IA),K=1,NWAVE)
  110   CONTINUE
        AMASS(IA)=AMASS(IA)*1822.89
  120 CONTINUE

*::::::::::::::::::::  initial wavefunction  ::::::::::::::::::::::::::
      OPEN(UNIT=14,FILE='ELCIN',STATUS='OLD',FORM='UNFORMATTED')
      READ(14,ERR=9140,END=9141) ICUBE,NFFT0,NWAVE0,UNIT0
      IF(ICUBE.NE.JCUBE) THEN
        IERR=110
        GO TO 9000
      ENDIF
      IF(NFFT0.NE.NFFT) THEN
        IERR=111
        GO TO 9000
      ENDIF
      IF(NWAVE0.NE.NWAVE) THEN
        IERR=112
        GO TO 9000
      ENDIF
      IF(UNIT0.NE.UNIT) THEN
        IERR=113
        GO TO 9000
      ENDIF
      READ(14,ERR=9140,END=9141) ISPIN0,NE
      IF(ISPIN0.NE.ISPIN) THEN
        IERR=114
        GO TO 9000
      ENDIF
      IF(ISPIN.EQ.1) NE(2)=0
      IF((NE(1).GT.NELC1).OR.(NE(2).GT.NELC2)) THEN
        WRITE(10,*) 'NE=',NE
        IERR=115
        GO TO 9000
      ENDIF
      N1(1)=1
      N2(1)=NE(1)
      N1(2)=NE(1)+1
      N2(2)=NE(1)+NE(2)
      DO 130 I=1,N2(ISPIN)
        READ(14,ERR=9140,END=9141) (C2(K,I),K=1,NWAVE)
  130 CONTINUE

      CALL SETUP(NI,ECUT,IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+200
        GO TO 9000
      ENDIF
*:::::::::::::::::  geometrical center of the cluster  ::::::::::::::::
      CX=dsum(NI,R2(1,1),1)/NI
      CY=dsum(NI,R2(1,2),1)/NI
      CZ=dsum(NI,R2(1,3),1)/NI

*:::::::::::::::::::::::::  center of mass  :::::::::::::::::::::::::::
      GX=0.0d0
      GY=0.0d0
      GZ=0.0d0
      AM=0.0d0
      DO 140 I=1,NI
        GX=GX+AMASS(KATM(I))*R2(I,1)
        GY=GY+AMASS(KATM(I))*R2(I,2)
        GZ=GZ+AMASS(KATM(I))*R2(I,3)
        AM=AM+AMASS(KATM(I))
  140 CONTINUE
      GX=GX/AM
      GY=GY/AM
      GZ=GZ/AM

*:::::::::::::::::::::  summary of input data  ::::::::::::::::::::::::
      WRITE(10,1110)
      WRITE(10,1115)

      IF(MOVE) THEN
        WRITE(10,1120) 'yes'
      ELSE
        WRITE(10,1120) 'no'
      ENDIF
      WRITE(10,1130) ISPIN
      WRITE(10,1140)
*     WRITE(10,1150)(I,ATOM(I),AMASS(I)/1822.89d0,ZV(I),LMAX(I),
*    &               I=1,NKATM)
      do ia = 1,nkatm
        write(10,1150) ia,atom(ia),amass(ia)/1822.89d0,
     >                 zv(ia),LMAX(ia)
        write(10,1151) (rc(i,ia),i=0,lmax(ia))
      end do
      WRITE(10,1160)
      WRITE(10,1170) (ATOM(K),NATOM(K),K=1,NKATM)
      WRITE(10,1180)
      WRITE(10,1190) (I,ATOM(KATM(I)),(R2(I,K),K=1,3),I=1,NI)
      WRITE(10,1200) CX,CY,CZ
      WRITE(10,1210) GX,GY,GZ
      WRITE(10,1220) NE(1),NE(ISPIN)
      WRITE(10,1230)
      IF(ICUBE.EQ.1) THEN
        WRITE(10,1240) 'simple cubic',UNIT,OMEGA
      ELSE IF(ICUBE.EQ.2) THEN
        WRITE(10,1240) 'face-centered cubic',UNIT,OMEGA
      ELSE
        WRITE(10,1240) 'body-centered cubic',UNIT,OMEGA
      ENDIF
      WRITE(10,1250) ECUT,NFFT,NFFT,NFFT,NIDA+2*NIDB
      WRITE(10,1260) RCUT,NCUT
      WRITE(10,1270)
      WRITE(10,1280) DT,FMASS
      WRITE(10,1290) TOLE,TOLC,TOLR
      WRITE(10,1300)
      WRITE(10,1305)
      CALL flush(10)
      CALL current_second(CPU2)
      ICOUNT=0
      E(1)=0

      IF(MOVE) OPEN(UNIT=30,FILE='MOTION',FORM='FORMATTED')

*:::::::::::::::::::::  begin iteration  ::::::::::::::::::::::::::::::
      CALL MESSAGE(2)
 1    CONTINUE
      ICOUNT=ICOUNT+1
      CALL CPSD(MOVE,DT,ITEST,E,DELTAE,DELTAC,DELTAR,IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+300
        GO TO 9000
      ENDIF

*:::::::::  store obtained wavefunctions and ionic geometry  ::::::::::
c     OPEN(UNIT=25,FILE='ELCBAK',FORM='UNFORMATTED')
c     OPEN(UNIT=26,FILE='IONBAK',FORM='FORMATTED')
c     WRITE(25) ICUBE,NFFT,NWAVE,UNIT
c     WRITE(25) ISPIN,NE
c     DO 180 I=1,N2(ISPIN)
c 180 WRITE(25) (C1(K,I),K=1,NWAVE)
c     WRITE(26,1370) NKATM
c     WRITE(26,1380)(FATOM(I),LMAX(I),I=1,NKATM)
c     WRITE(26,1390) NI
c     WRITE(26,1400) (KATM(J),(R2(J,I),I=1,3),(R1(J,I),I=1,3),J=1,NI)
c     CLOSE(25)
c     CLOSE(26)

      WRITE(10,1310) ICOUNT*ITEST,E(1),DELTAE,DELTAC,DELTAR
      CALL flush(10)
      IF(MOVE) THEN
        WRITE(30,'(I7)') ICOUNT
        WRITE(30,'(3F11.5)') (R2(I,1),R2(I,2),R2(I,3),I=1,NI)
      ENDIF
      IF(DELTAE.GT.0.) THEN
        WRITE(10,*) ' *** ENERGY GOING UP.  iteration terminated.'
        GO TO 2
      ENDIF
      DELTAE=dabs(DELTAE)
      IF(DELTAE.LT.TOLE.AND.DELTAC.LT.TOLC.AND.DELTAR.LT.TOLR) THEN
        WRITE(10,*) ' *** tolerences ok.   iteration terminated.'
        GO TO 2
      ENDIF
      IF(ICOUNT.LT.ISTOP) GO TO 1
      WRITE(10,*) '*** arrived at the MAXIMUM iteration.   terminated.'
*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

    2 CONTINUE
      CALL MESSAGE(3)
      CALL current_second(CPU3)

*::::::::::::::  check total number of electrons  :::::::::::::::::::::
      OPEN(UNIT=17,FILE='CHECK',STATUS='NEW',FORM='FORMATTED')
      DO 200 MS=1,ISPIN
        EN(MS)=dsum(2*NFFT3D,DN(1,MS),1)*OMEGA/dble(NFFT**3)
  200 CONTINUE
      WRITE(17,1320) (EN(MS),MS=1,ISPIN)

*::::::::: comparison between hamiltonian and lambda matrix :::::::::::
      WRITE(17,1330)
      DO 210 MS=1,ISPIN
        DO 210 I=1,NE(MS)
          DO 210 J=I,NE(MS)
             WRITE(17,1340) MS,I,J,HML(I,J,MS),LMD(J,I,MS),
     &                      HML(I,J,MS)-LMD(J,I,MS)
  210 CONTINUE

*:::::::::::::::::::::  check orthnormality  ::::::::::::::::::::::::::
      WRITE(17,1350)
      DO 220 MS=1,ISPIN
        DO 220 I=N1(MS),N2(MS)
          II=I-N1(MS)+1
          DO 220 J=I,N2(MS)
            JJ=J-N1(MS)+1
            W=GCDOTC(NIDA,NIDB,C1(1,I),C1(1,J))
            WRITE(17,1360) MS,II,JJ,W
  220 CONTINUE

*::::::::::::::::  diagnalization of hamiltonian matrix  ::::::::::::::
      CALL dcopy(ISPIN*NELC1,0.0d0,0,EIG,1)
      DO 225 MS=1,ISPIN
        CALL EIGEN(NELC1,NE(MS),HML(1,1,MS),EIG(1,MS),XTMP,IERR)
        IF(IERR.NE.0) THEN
          IERR=IERR+400
          GO TO 9000
        ENDIF
  225 CONTINUE
      CALL dcopy(2*N2(ISPIN)*NWAVE,0.0d0,0,C1,1)
      DO 230 MS=1,ISPIN
        DO 230 J=N1(MS),N2(MS)
          JJ=J-N1(MS)+1
          DO 230 I=N1(MS),N2(MS)
            II=I-N1(MS)+1
            CALL daxpy(2*NWAVE,HML(II,JJ,MS),C2(1,I),1,C1(1,J),1)
  230 CONTINUE

*:::::::::  store obtained wavefunctions and ionic geometry  ::::::::::
      OPEN(UNIT=15,FILE='ELCOUT', STATUS='NEW',FORM='UNFORMATTED')
      OPEN(UNIT=16,FILE='IONOUT', STATUS='NEW',FORM='FORMATTED')
      OPEN(UNIT=18,FILE='DNSOUT', STATUS='NEW',FORM='UNFORMATTED')
      WRITE(15) ICUBE,NFFT,NWAVE,UNIT
      WRITE(15) ISPIN,NE
      DO 240 I=1,N2(ISPIN)
        WRITE(15) (C1(K,I),K=1,NWAVE)
  240 CONTINUE
      DO 250 I=1,NI
        R1(I,1)=0.0d0
        R1(I,2)=0.0d0
        R1(I,3)=0.0d0
  250 CONTINUE
      WRITE(16,1370) NKATM
      WRITE(16,1380)(FATOM(I),LMAX(I),I=1,NKATM)
      WRITE(16,1390) NI
      WRITE(16,1400) (KATM(J),(R2(J,I),I=1,3),(R1(J,I),I=1,3),J=1,NI)

      WRITE(18) ICUBE,NFFT,UNIT
      WRITE(18) ISPIN
      DO 252 I=1,ISPIN
        WRITE(18) (DN(K,I),K=1,2*NFFT3D)
  252 CONTINUE
      close(15)
      close(16)
      close(18)

*:::::::::::::::::  geometrical center of the cluster  ::::::::::::::::
      CX=dsum(NI,R2(1,1),1)/NI
      CY=dsum(NI,R2(1,2),1)/NI
      CZ=dsum(NI,R2(1,3),1)/NI

*:::::::::::::::::::::::::  center of mass  :::::::::::::::::::::::::::
      GX=0.0d0
      GY=0.0d0
      GZ=0.0d0
      AM=0.0d0
      DO 260 I=1,NI
        GX=GX+AMASS(KATM(I))*R2(I,1)
        GY=GY+AMASS(KATM(I))*R2(I,2)
        GZ=GZ+AMASS(KATM(I))*R2(I,3)
        AM=AM+AMASS(KATM(I))
  260 CONTINUE
      GX=GX/AM
      GY=GY/AM
      GZ=GZ/AM

*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      WRITE(10,1300)
      WRITE(10,1410)
      WRITE(10,1420)
      WRITE(10,1190)(I,ATOM(KATM(I)),(R2(I,K),K=1,3),I=1,NI)
      WRITE(10,1200) CX,CY,CZ
      WRITE(10,1210) GX,GY,GZ
      WRITE(10,1430) E(1),E(1)/NI
      WRITE(10,1440) E(2),E(2)/N2(ISPIN)
      WRITE(10,1450) E(3),E(3)/N2(ISPIN)
      WRITE(10,1460) E(4),E(4)/N2(ISPIN)
      WRITE(10,1470) E(5),E(5)/NI
      WRITE(10,1480) E(6),E(6)/N2(ispin)
      WRITE(10,1490) E(7),E(7)/N2(ispin)
      WRITE(10,1495) E(8),E(8)/N2(ispin)
      WRITE(10,1496) E(9),E(9)/N2(ispin)
      WRITE(10,1497) E(10),E(10)/N2(ispin)
      virial = (E(10)+E(9)+E(8)+E(7))/E(6)
      WRITE(10,1498) virial
      WRITE(10,1500)
      NN=NE(1)-NE(2)
      EV=27.2116
      DO 270 I=1,NN
        WRITE(10,1510) EIG(I,1),EIG(I,1)*EV
  270 CONTINUE
      DO 280 I=1,NE(2)
        WRITE(10,1510) EIG(I+NN,1),EIG(I+NN,1)*EV,EIG(I,2),EIG(I,2)*EV
  280 CONTINUE

      CALL current_second(CPU4)

*:::::::::::::::::::   report consumed cputime  :::::::::::::::::::::::
      T1=CPU2-CPU1
      T2=CPU3-CPU2
      T3=CPU4-CPU3
      T4=CPU4-CPU1
      AV=T2/(ICOUNT*ITEST)
      WRITE(10,*)
      WRITE(10,*) '-----------------'
      WRITE(10,*) 'cputime in seconds'
      WRITE(10,*) 'prologue    : ',T1
      WRITE(10,*) 'main loop   : ',T2
      WRITE(10,*) 'epilogue    : ',T3
      WRITE(10,*) 'total       : ',T4
      WRITE(10,*) 'cputime/step: ',AV
      CALL MESSAGE(4)
      STOP
*::::::::::::::::::::  I/O error codes  ::::::::::::::::::::::::::::::
 9110 IERR=10
      GO TO 9000
 9111 IERR=11
      GO TO 9000
 9120 IERR=12
      GO TO 9000
 9121 IERR=13
      GO TO 9000
 9130 IERR=14
      GO TO 9000
 9131 IERR=15
      GO TO 9000
 9140 IERR=16
      GO TO 9000
 9141 IERR=17
      GO TO 9000

*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*         [ steepest descent minimization ]        *')
 1035 FORMAT(10x,'*         [ Double precision on ibm590    ]        *')
 1040 FORMAT(10X,'*            version #3.00   09/10/89              *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1115 FORMAT(/' options:')
 1120 FORMAT(5X,' ionic motion   = ',A)
 1130 FORMAT(5X,' electron spin  = ',I1)
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A2,'  mass no.:',F6.1,'   core charge:',F4.0,
     &       '  LMAX=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I2,'  spin down=',I2)
 1230 FORMAT(/' SUPERCELL:')
 1240 FORMAT(5X,' lattice: ',A,'   size=',F7.2,'  volume=',F10.1)
 1250 FORMAT(5X,' cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves)')
 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictacious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1300 FORMAT(//)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/ion)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l energy     :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Ratio <V>/<T>:',E19.10)
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 WRITE(10,9010) IERR
      CALL MESSAGE(5)
      STOP
      END
