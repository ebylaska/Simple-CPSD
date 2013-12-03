      SUBROUTINE CPSD(MOVE,DT,MXTIME,E,DELTAE,DELTAC,DELTAR,IERR)

      IMPLICIT double precision (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      LOGICAL MOVE
      double precision LMD
      COMPLEX*16 S,DNG,C2,C1,CPSI,VALL,VTMP,CFFT,CTMP
      COMPLEX*16 EX1,EX2,EX3,EXI
      COMPLEX*16 CW1,CW2,CW3

      INCLUDE 'aimd.h'

*:::::::::  secondary parameters that determine array size ::::::::::::
*     NEMAX  : maximum number of electrons
*     NFFT3D : size of array for fast fourier transform
*     NSHL3D : number of lattice to be summed in the Ewald summation
      PARAMETER (NEMAX=NELC1+NELC2)
      PARAMETER (NFFT3D=(NFFT/2+1)*NFFT*NFFT,N2FT3D=2*NFFT3D)
      PARAMETER (NWAVE2=2*NWAVE)
      PARAMETER (NSHL3D=(2*NSHELL+1)**3)


      PARAMETER (M1=NFFT,M2=NFFT,M3=NFFT,M0=M1+2)
      PARAMETER (M1H=M1/2,M2H=M2/2,M3H=M3/2,M0H=M0/2)
      PARAMETER (INC2R=M0,INC3R=M0*M2,INC2C=M0H,INC3C=M0H*M2)
      PARAMETER (IFWD=1,IBWD=-1)

*::::::::::::::::::  constants for MFFT routines  :::::::::::::::::::::
*     parameter (IOPT=1,IORD=1,ISET=0)
*     real      mw1,mw2,mw3,iwork
*     dimension mw1(6*M1+14),mw2(4*M2*(M0/2+1)+14),mw3(4*M3+14)
*     dimension iwork(M1)
*
*:::::::::::::::::::::: Constants for Essl ffts :::::::::::::::::::::::
      parameter (naux=100000)
      double precision aux(naux)


*:::::::::::  expansion coefficient of the error function  ::::::::::::
      PARAMETER (CERFC=1.128379167d0)
      PARAMETER (B1=0.0705230784d0,B2=0.0422820123d0,B3=0.0092705272d0)
      PARAMETER (B4=0.0001520143d0,B5=0.0002765672d0,B6=0.0000430638d0)

*:::::::  iteration limit and tolerence for non-liner equations  ::::::
      PARAMETER (ITRLMD=20, CONVG=1.d-10)

      DIMENSION KATM(NATMX),VL(NWAVE,NKATMX),VNL(NWAVE,9,NKATMX)
      DIMENSION ZV(NKATMX),LMMAX(NKATMX),AMASS(NKATMX),VNLNRM(9,NKATMX)
      DIMENSION INDX(NWAVE)
      DIMENSION VC(NWAVE),VG(NWAVE),TG(NWAVE),G(NWAVE,3)
      DIMENSION UNITA(3,3),UNITG(3,3),RCELL(NSHL3D,3),FT(NSHL3D,3)
      DIMENSION NE(2),N1(2),N2(2),E(10)
      DIMENSION R2(NATMX,3),R1(NATMX,3),FI(NATMX,3),DTI(NATMX)
      DIMENSION DN(N2FT3D,2),XCP(N2FT3D,2),XCE(N2FT3D)
      DIMENSION LMD(NELC1,NELC1,2),HML(NELC1,NELC1,2)
      DIMENSION A22(NELC1,NELC1),A21(NELC1,NELC1),A12(NELC1,NELC1)
      DIMENSION A11(NELC1,NELC1),AA0(NELC1,NELC1),AA1(NELC1,NELC1)
      DIMENSION S(NWAVE),DNG(NWAVE)
      DIMENSION SPSI(N2FT3D,NEMAX),CPSI(NFFT3D,NEMAX)
      DIMENSION SFFT(N2FT3D),CFFT(NFFT3D),XTMP(N2FT3D)
      DIMENSION VTMP(NWAVE),CTMP(NFFT3D),STMP(N2FT3D)
      DIMENSION EXI(NWAVE),VALL(NWAVE)
      DIMENSION C2(NWAVE,NEMAX),C1(NWAVE,NEMAX)
      DIMENSION EX1(0:M1H,NATMX),EX2(0:M2-1,NATMX),EX3(0:M3-1,NATMX)


      EQUIVALENCE (SPSI,CPSI),(SFFT,CFFT),(STMP,CTMP)

      COMMON            c2,c1,lmd,hml,fmass,ne
      COMMON / ION    / r2,r1,amass,ni
      COMMON / DENSTY / dn,xcp,xce,xtmp,ispin
      COMMON / FOURIE / tg,vc,vg,g
      COMMON / LATTIC / unita,unitg,unit,omega,icube
      COMMON / EWALD  / rcell,cewald,rcut,ncut
      COMMON / PSEUDO / vl,vnl,vnlnrm,zv,nkatm,katm,lmmax
      COMMON / INDEX2  / indx,nida,nidb
      COMMON / WORKSP / stmp



*:::::::::::::::::::::::::::  CONSTANTS  :::::::::::::::::::::::::::::::
      PI=4.0d0*DATAN(1.0d0)
      N1(1)=1
      N2(1)=NE(1)
      N1(2)=NE(1)+1
      N2(2)=NE(1)+NE(2)
      SCAL1=1.0d0/dble(M1*M2*M3)
      SCAL2=1.0d0/OMEGA
      DV=OMEGA*SCAL1
      DTE=DT/dSQRT(FMASS)
      IF(MOVE) THEN
        DO 110 I=1,NI
          DTI(I)=DT/dSQRT(AMASS(KATM(I)))
  110   CONTINUE
      ENDIF


*::::::::::::::::  initialization of MFFT routines  :::::::::::::::::::
*     CALL R3FFT(CFFT,M0,M1,M2,M3,MW1,MW2,MW3,IOPT,ISET,IORD,IWORK,IERR)
*     IF(IERR.NE.0) RETURN

*                       =====================
*::::::::::::::::::::   |  start main loop  |   ::::::::::::::::::::::::
*                       =====================
      ITIME=0
  200 CONTINUE

*::::::::::::::::::::::  increment time  ::::::::::::::::::::::::::::::
      ITIME=ITIME+1
      CALL dcopy(N2(ISPIN)*NWAVE2,C2,1,C1,1)
      CALL dcopy(NI,R2(1,1),1,R1(1,1),1)
      CALL dcopy(NI,R2(1,2),1,R1(1,2),1)
      CALL dcopy(NI,R2(1,3),1,R1(1,3),1)

*::::::::::::::  wavefunctions in the coordination space  ::::::::::::::

      call dcopy(N2(ISPIN)*N2FT3D,0.0d0,0,SPSI,1)
      DO 210 N=N1(1),N2(ISPIN)
        DO 205 K=1,NWAVE
          CPSI(INDX(K),N)=C1(K,N)
  205   CONTINUE

*:::::::::::::::::::::: Cray :::::::::::::::::::::::::::::::::::::
*       CALL R3FFT(CPSI(1,N),M0,M1,M2,M3,MW1,MW2,MW3,IOPT,IBWD,IORD,
*    &             IWORK,IERR)
*::::::::::::::::::::::: Essl - complex to real fft :::::::::::::::
*       call dcrft3(cpsi(1,n), inc2c, inc3c,
*    >              spsi(1,n), inc2r, inc3r,
*    >	            nfft, nfft,nfft,ibwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  complex to real fft :::::::::::::::
        call cr_fft3b(cpsi(1,n),
     >                 inc2c,inc3c,inc2r,inc3r,
     >	               nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  210 CONTINUE

*:::::::::::::::::::::::  electron spin density  ::::::::::::::::::::::
      call dcopy(ISPIN*N2FT3D,0.0d0,0,DN,1)
      do MS=1,ISPIN
        do N=N1(MS),N2(MS)
          do K=1,N2FT3D
            DN(K,MS)=DN(K,MS)+SCAL2*SPSI(K,N)**2
          end do
        end do
        call dcopy(NFFT*NFFT,0.0d0,0,dn(M1+1,ms),INC2R)
        call dcopy(NFFT*NFFT,0.0d0,0,dn(M1+2,ms),INC2R)
        sw1 = dsum(n2ft3d,dn(1,ms),1)*OMEGA/dble(NFFT**3)
        write(*,*) "DN:",sw1
      end do
      

*:::::::::::::  fourier transform of the electron density  ::::::::::::
      DO 230 K=1,N2FT3D
        SFFT(K)=(DN(K,1)+DN(K,ISPIN))*SCAL1
  230 CONTINUE
      CALL dcopy(NFFT*NFFT,0.0d0,0,SFFT(M1+1),INC2R)
      CALL dcopy(NFFT*NFFT,0.0d0,0,SFFT(M1+2),INC2R)

*:::::::::::::::::::::: Cray :::::::::::::::::::::::::::::::::::::
*     CALL R3FFT(SFFT,M0,M1,M2,M3,MW1,MW2,MW3,IOPT,IFWD,IORD,IWORK,IERR)
*
*::::::::::::::::::: Essl - Real to complex ::::::::::::::::::::::::::
*     call drcft3(sfft, inc2r, inc3r,
*    >            cfft, inc2c, inc3c,
*    >            nfft,nfft,nfft,ifwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  real to complex fft :::::::::::::::
        call rc_fft3f(cfft,
     >                 inc2c,inc3c,inc2r,inc3r,
     >	               nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      DO 235 K=1,NWAVE
        DNG(K)=CFFT(INDX(K))
  235 CONTINUE

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
      DO 270 I=1,NI
        SW1=UNITG(1,1)*R1(I,1)+UNITG(2,1)*R1(I,2)+UNITG(3,1)*R1(I,3)+PI
        SW2=UNITG(1,2)*R1(I,1)+UNITG(2,2)*R1(I,2)+UNITG(3,2)*R1(I,3)+PI
        SW3=UNITG(1,3)*R1(I,1)+UNITG(2,3)*R1(I,2)+UNITG(3,3)*R1(I,3)+PI
        CW1=dcmplx(dcos(SW1),-dsin(SW1))
        CW2=dcmplx(dcos(SW2),-dsin(SW2))
        CW3=dcmplx(dcos(SW3),-dsin(SW3))
        EX1(0,I)=dcmplx(1.0d0,0.0d0)
        EX2(0,I)=dcmplx(1.0d0,0.0d0)
        EX3(0,I)=dcmplx(1.0d0,0.0d0)
        DO 240 K=1,M1H
          EX1(K,I)=EX1(K-1,I)*CW1
  240   CONTINUE
        DO 250 K=1,M2H
          EX2(K,I)=EX2(K-1,I)*CW2
          EX2(M2-K,I)=dconjg(EX2(K,I))
  250   CONTINUE
        DO 260 K=1,M3H
          EX3(K,I)=EX3(K-1,I)*CW3
          EX3(M3-K,I)=dconjg(EX3(K,I))
  260   CONTINUE
c       EX1(M1H,I)=dcmplx(2.0d0*dble(EX1(M1H,I)),0.0d0)
c       EX2(M2H,I)=dcmplx(2.0d0*dble(EX2(M2H,I)),0.0d0)
c       EX3(M3H,I)=dcmplx(2.0d0*dble(EX3(M3H,I)),0.0d0)
        EX1(M1H,I)=0.0d0
        EX2(M2H,I)=0.0d0
        EX3(M3H,I)=0.0d0
  270 CONTINUE

*:::::::::::::::::::::::::  kinetic energy  :::::::::::::::::::::::::::
      DO 280 N=N1(1),N2(ISPIN)
        DO 280 K=1,NWAVE
          C2(K,N)=-TG(K)*C1(K,N)
  280 CONTINUE


*::::::::::::::::::  start potential energy part  :::::::::::::::::::::

      CALL dcopy(NWAVE2,0.0d0,0,S,1)
      CALL dcopy(NWAVE2,0.0d0,0,VALL,1)

      DO 380 I=1,NI
        IA=KATM(I)

*::::::::::::::::  structure factor and local pseudopotential  ::::::::
        DO 290 K3=0,M3-1
          DO 290 K2=0,M2-1
            DO 290 K1=0,M1H
            CFFT(INC3C*K3+INC2C*K2+K1+1)=EX1(K1,I)*EX2(K2,I)*EX3(K3,I)
  290   CONTINUE
        DO 295 K=1,NWAVE
          EXI(K)=CFFT(INDX(K))
  295   CONTINUE
        CALL DAXPY(NWAVE2,ZV(IA),EXI,1,S,1)
        DO 300 K=1,NWAVE
          VTMP(K)=VL(K,IA)*EXI(K)
          VALL(K)=VALL(K)+VTMP(K)
  300   CONTINUE
        IF(MOVE) THEN
          DO 310 K=1,NWAVE
            XTMP(K)=+dimag(DNG(K))*dble(VTMP(K))
     &              -dble(DNG(K))*dimag(VTMP(K))
  310     CONTINUE
          FI(I,1)=GSDOT(NIDA,NIDB,G(1,1),XTMP)
          FI(I,2)=GSDOT(NIDA,NIDB,G(1,2),XTMP)
          FI(I,3)=GSDOT(NIDA,NIDB,G(1,3),XTMP)
       ENDIF
       

*::::::::::::::::::::  non-local pseudopotential  :::::::::::::::::::::
        DO 370 L=1,LMMAX(IA)
          IF(L.EQ.1 .OR. L.GE.5) THEN
            DO 330 K=1,NWAVE
              VTMP(K)=VNL(K,L,IA)*EXI(K)
  330       CONTINUE
          ELSE
            DO 335 K=1,NWAVE
              VTMP(K)=VNL(K,L,IA)*dcmplx(-dimag(EXI(K)),
     >                                     dble(EXI(K)))
  335       CONTINUE
          ENDIF
          DO 360 N=N1(1),N2(ISPIN)
            DO 340 K=1,NWAVE
              CTMP(K)=C1(K,N)*dconjg(VTMP(K))
  340       CONTINUE
            SW1=GCSUM(NIDA,NIDB,CTMP)/OMEGA/VNLNRM(L,IA)
            CALL daxpy(NWAVE2,-SW1,VTMP,1,C2(1,N),1)
            IF(MOVE) THEN
              IF(ISPIN.EQ.1) SW1=2.0d0*SW1
              CALL dcopy(NWAVE,STMP(2),2,XTMP,1)
              FI(I,1)=FI(I,1)+2.0d0*SW1*GSDOT(NIDA,NIDB,G(1,1),XTMP)
              FI(I,2)=FI(I,2)+2.0d0*SW1*GSDOT(NIDA,NIDB,G(1,2),XTMP)
              FI(I,3)=FI(I,3)+2.0d0*SW1*GSDOT(NIDA,NIDB,G(1,3),XTMP)
            ENDIF
  360     CONTINUE
  370   CONTINUE
        
  380 CONTINUE
      

*::::::::::::::::::::::::  Ewald summation  :::::::::::::::::::::::::::
      IF(MOVE) THEN
        DO 420 I=1,NI
          IA=KATM(I)
          DO 390 K3=0,M3-1
            DO 390 K2=0,M2-1
              DO 390 K1=0,M1H
                CFFT(INC3C*K3+INC2C*K2+K1+1)=EX1(K1,I)*EX2(K2,I)
     &                                      *EX3(K3,I)
  390     CONTINUE
          DO 395 K=1,NWAVE
            EXI(K)=CFFT(INDX(K))
  395     CONTINUE
          DO 400 K=1,NWAVE
              XTMP(K)=(dble(EXI(K))*dimag(S(K))
     &               -dimag(EXI(K))*dble(S(K)))*VG(K)
  400     CONTINUE
          FI(I,1)=FI(I,1)+GSDOT(NIDA,NIDB,G(1,1),XTMP)*ZV(IA)*SCAL2
          FI(I,2)=FI(I,2)+GSDOT(NIDA,NIDB,G(1,2),XTMP)*ZV(IA)*SCAL2
          FI(I,3)=FI(I,3)+GSDOT(NIDA,NIDB,G(1,3),XTMP)*ZV(IA)*SCAL2
  420 CONTINUE

        DO 440 I=1,NI
          DO 440 J=I+1,NI
            DX=R1(I,1)-R1(J,1)
            DY=R1(I,2)-R1(J,2)
            DZ=R1(I,3)-R1(J,3)
            ZZ=ZV(KATM(I))*ZV(KATM(J))
            DO 430 L=1,NSHL3D
              X=RCELL(L,1)+DX
              Y=RCELL(L,2)+DY
              Z=RCELL(L,3)+DZ
              R=dsqrt(X*X+Y*Y+Z*Z)
              W=R/RCUT
              ERFC=(1.0d0+W*(B1+W*(B2+W*(B3+W*(B4
     >                   +W*(B5+W*B6))))))**4
              ERFC=1.0d0/ERFC**4
              F=ZZ*(ERFC+CERFC*W*DEXP(-W*W))/R**3
              FT(L,1)=X*F
              FT(L,2)=Y*F
              FT(L,3)=Z*F
  430       CONTINUE
            SW1=dsum(NSHL3D,FT(1,1),1)
            SW2=dsum(NSHL3D,FT(1,2),1)
            SW3=dsum(NSHL3D,FT(1,3),1)
            FI(I,1)=FI(I,1)+SW1
            FI(I,2)=FI(I,2)+SW2
            FI(I,3)=FI(I,3)+SW3
            FI(J,1)=FI(J,1)-SW1
            FI(J,2)=FI(J,2)-SW2
            FI(J,3)=FI(J,3)-SW3
  440   CONTINUE
      ENDIF
      


*:::::  local pseudo-, hartree-, exchange-correlation potentials  :::::

      CALL VXC

      call dcopy(N2FT3D,0.0d0,0,SFFT,1)
      DO 450 K=1,NWAVE
        CFFT(INDX(K))=VALL(K)*SCAL2+VC(K)*DNG(K)
  450 CONTINUE

*:::::::::::::::::::::: Cray :::::::::::::::::::::::::::::::::::::
*     CALL R3FFT(CFFT,M0,M1,M2,M3,MW1,MW2,MW3,IOPT,IBWD,IORD,IWORK,IERR)
*     IF(IERR.NE.0) RETURN
*     CALL dcopy(N2FT3D,CFFT,1,XTMP,1)
*
*::::::::::::::::::::::: Essl - complex to real fft :::::::::::::::
*       call dcrft3(cfft, inc2c, inc3c,
*    >              xtmp, inc2r, inc3r,
*    >	            nfft, nfft,nfft,ibwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  complex to real fft :::::::::::::::
        call cr_fft3b(cfft,
     >                 inc2c,inc3c,inc2r,inc3r,
     >	               nfft,nfft,nfft,aux,naux)
      call dcopy(n2ft3d,cfft,1,xtmp,1)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



      DO 490 MS=1,ISPIN
        DO 460 K=1,N2FT3D
          STMP(K)=XTMP(K)+XCP(K,MS)
  460   CONTINUE
        DO 480 N=N1(MS),N2(MS)
          DO 470 K=1,N2FT3D
            SFFT(K)=STMP(K)*SPSI(K,N)
  470     CONTINUE


*:::::::::::::::::::::: Cray :::::::::::::::::::::::::::::::::::::
*         CALL R3FFT(SFFT,M0,M1,M2,M3,MW1,MW2,MW3,IOPT,IFWD,IORD,
*    &               IWORK,IERR)
*
*::::::::::::::::::: Essl - Real to complex ::::::::::::::::::::::::::
*     call drcft3(sfft, inc2r, inc3r,
*    >            cfft, inc2c, inc3c,
*    >            nfft,nfft,nfft,ifwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  real to complex fft :::::::::::::::
        call rc_fft3f(cfft,
     >                 inc2c,inc3c,inc2r,inc3r,
     >	               nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

         do k=1,nwave
            cpsi(k,n)=cfft(indx(k))
         end do

*     store vector H|PSI> in the array CPSI
          DO 475 K=1,NWAVE
            CPSI(K,N)=-SCAL1*CFFT(INDX(K))
            CPSI(K,N)=CPSI(K,N)+C2(K,N)
  475     CONTINUE

  480   CONTINUE
  490 CONTINUE


*::::::::::::::::::::  steepest descent equations  ::::::::::::::::::::
      DO 560 N=1,N2(ISPIN)
*       CALL dzaxpy(NWAVE2,DTE,CPSI(1,N),1,C1(1,N),1,C2(1,N),1)
        call dcopy(nwave2,c1(1,n),1,c2(1,n),1)
        call daxpy(NWAVE2,DTE,CPSI(1,N),1,C2(1,N),1)
  560 CONTINUE

      IF(MOVE) THEN
        DO 570 I=1,NI
          R2(I,1)=R1(I,1)+DTI(I)*FI(I,1)
          R2(I,2)=R1(I,2)+DTI(I)*FI(I,2)
          R2(I,3)=R1(I,3)+DTI(I)*FI(I,3)
  570   CONTINUE
      ENDIF

*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      DO 640 MS=1,ISPIN
        IF(NE(MS).LE.0) GO TO 640
        DO 610 II=N1(MS),N2(MS)
          I=II-N1(MS)+1
          A22(I,I)=(1.0d0-GCDOTC(NIDA,NIDB,C2(1,II),C2(1,II)))
     >                *0.5d0/DTE
          A21(I,I)=(1.0d0-GCDOTC(NIDA,NIDB,C2(1,II),C1(1,II)))
     >                *0.5d0
          A12(I,I)=A21(I,I)
          A11(I,I)=-GCDOTC(NIDA,NIDB,C1(1,II),C1(1,II))
     >                 *0.5d0*DTE

          DO 610 JJ=II+1,N2(MS)
            J=JJ-N1(MS)+1
            A22(I,J)=-GCDOTC(NIDA,NIDB,C2(1,II),C2(1,JJ))
     >                   *0.5d0/DTE
            A21(I,J)=-GCDOTC(NIDA,NIDB,C2(1,II),C1(1,JJ))*0.5d0
            A12(I,J)=-GCDOTC(NIDA,NIDB,C1(1,II),C2(1,JJ))*0.5d0
            A11(I,J)=-GCDOTC(NIDA,NIDB,C1(1,II),C1(1,JJ))
     >                 *0.5d0*DTE
            A22(J,I)=A22(I,J)
            A21(J,I)=A12(I,J)
            A12(J,I)=A21(I,J)
            A11(J,I)=A11(I,J)
 610    CONTINUE
        CALL dcopy(NELC1*NE(MS),A22,1,AA0,1)
        DO 620 IT=1,ITRLMD
          CALL dcopy(NELC1*NE(MS),A22,1,AA1,1)
          CALL DMMUL(NELC1,NE(MS), A21, AA0,SFFT)
          CALL DMMUL(NELC1,NE(MS), AA0, A12,STMP)
          CALL DMADD(NELC1,NE(MS),SFFT,STMP,SFFT)
          CALL DMADD(NELC1,NE(MS),SFFT, AA1, AA1)
          CALL DMMUL(NELC1,NE(MS), A11, AA0,SFFT)
          CALL DMMUL(NELC1,NE(MS), AA0,SFFT,STMP)
          CALL DMADD(NELC1,NE(MS),STMP, AA1, AA1)
          CALL DMSUB(NELC1,NE(MS), AA1, AA0,SFFT)
          ADIFF=SFFT(idamax(NELC1*NE(MS),SFFT,1))
          IF(ADIFF.LT.CONVG) GO TO 630
          CALL dcopy(NELC1*NE(MS),AA1,1,AA0,1)
  620   CONTINUE
        IERR=10
        WRITE(6,*) 'IERR=10',ADIFF
C       RETURN
  630   CONTINUE
        CALL dcopy(NELC1*NE(MS),AA1,1,LMD(1,1,MS),1)
  640 CONTINUE

*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      DO 660 MS=1,ISPIN
        DO 650 II=N1(MS),N2(MS)
          I=II-N1(MS)+1
          DO 650 JJ=N1(MS),N2(MS)
            J=JJ-N1(MS)+1
            CALL daxpy(NWAVE2,DTE*LMD(I,J,MS),C1(1,JJ),1,C2(1,II),1)
  650   CONTINUE
  660 CONTINUE


      IF(ITIME.LT.MXTIME) GO TO 200
*::::::::::::::::::::::  end of main loop  ::::::::::::::::::::::::::::



*                    ========================
*::::::::::::::::::  total energy calculation  ::::::::::::::::::::::::
*                    ========================

*:::::::::::  hamiltonian matrix and orbital energy  ::::::::::::::::::
      DO 670 MS=1,ISPIN
        DO 670 II=N1(MS),N2(MS)
          I=II-N1(MS)+1
          HML(I,I,MS)=-GCDOTC(NIDA,NIDB,C1(1,II),CPSI(1,II))
          DO 670 JJ=II+1,N2(MS)
            J=JJ-N1(MS)+1
            HML(I,J,MS)=-GCDOTC(NIDA,NIDB,C1(1,II),CPSI(1,JJ))
            HML(J,I,MS)=HML(I,J,MS)
  670 CONTINUE
      EORBIT=0.0
      DO 675 MS=1,ISPIN
        DO 675 I=1,NE(MS)
          EORBIT=EORBIT+HML(I,I,MS)
  675 CONTINUE
      IF(ISPIN.EQ.1) EORBIT=EORBIT+EORBIT

*:::::::::::  ion-ion interaction energy (Ewald summation)  :::::::::::
      DO 700 K=1,NWAVE
        XTMP(K)=(dble(S(K))**2+dimag(S(K))**2)*VG(K)
  700 CONTINUE
      EION=GSSUM(NIDA,NIDB,XTMP)*0.5d0/OMEGA+CEWALD
      
      DO 705 I=1,NI
        DO 705 J=I+1,NI
          DX=R1(I,1)-R1(J,1)
          DY=R1(I,2)-R1(J,2)
          DZ=R1(I,3)-R1(J,3)
          ZZ=ZV(KATM(I))*ZV(KATM(J))
          DO 702 L=1,NSHL3D
            X=RCELL(L,1)+DX
            Y=RCELL(L,2)+DY
            Z=RCELL(L,3)+DZ
            R=dSQRT(X*X+Y*Y+Z*Z)
            W=R/RCUT
            ERFC=1.0d0/(1.0d0+W*(B1+W*(B2+W*(B3
     >                   +W*(B4+W*(B5+W*B6))))))**4
            FT(L,1)=ZZ*ERFC**4/R
  702     CONTINUE
          EION=EION+dsum(NSHL3D,FT(1,1),1)
  705 CONTINUE
      

*::::::::::::::::::::::  Hartree energy  ::::::::::::::::::::::::::::::
      DO 710 K=1,NWAVE
        XTMP(K)=(dble(DNG(K))**2+dimag(DNG(K))**2)*VC(K)
  710 CONTINUE
      EHARTR=0.5d0*GSSUM(NIDA,NIDB,XTMP)*OMEGA

*:::::::::::::::::::::  exchange-correlation energy  ::::::::::::::::::
      EXC=ddot(N2FT3D,DN(1,1),1,XCE,1)
      PXC=ddot(N2FT3D,DN(1,1),1,XCP(1,1),1)
      IF(ISPIN.EQ.1) THEN
        EXC=EXC+EXC
        PXC=PXC+PXC
      ELSE
        EXC=EXC+ddot(N2FT3D,DN(1,2),1,XCE,1)
        PXC=PXC+ddot(N2FT3D,DN(1,2),1,XCP(1,2),1)
      ENDIF
      EXC=EXC*DV
      PXC=PXC*DV

*     ***** Kohn-Sham kinetic energy ****
      eke = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do k=1,nwave
               cfft(k) = tg(k)*c1(k,n)
            end do
            eke = eke + gcdotc(nida,nidb,c1(1,n),cfft)
         end do
      end do
      if (ispin.eq.1) eke = 2.0d0*eke

*     **** Kohn-Sham V_local energy ****
      call dcopy(n2ft3d,0.0d0,0,cfft,1)
      do k=1,nwave
         cfft(indx(k)) = vall(k)*scal2
      end do
*::::::::::::::::  ESSL      -  complex to real fft :::::::::::::::
*     call dcrft3(cfft,inc2c,inc3c,
*    >            sfft,inc2r,inc3r,
*    >            nfft,nfft,nfft,ibwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  complex to real fft :::::::::::::::
        call cr_fft3b(cfft,
     >                 inc2c,inc3c,inc2r,inc3r,
     >	               nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      elocal = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call dcopy(n2ft3d,0.0d0,0,ctmp,1)
            do k=1,nwave
               ctmp(indx(k)) = c1(k,n)
            end do
*::::::::::::::::  ESSL      -  complex to real fft :::::::::::::::
*           call dcrft3(ctmp,inc2c,inc3c,
*    >                  stmp,inc2r,inc3r,
*    >                  nfft,nfft,nfft,ibwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  complex to real fft :::::::::::::::
            call cr_fft3b(ctmp,
     >                    inc2c,inc3c,inc2r,inc3r,
     >	                  nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            do k=1,n2ft3d
               stmp(k) = sfft(k)*stmp(k)
            end do
*::::::::::::::::  ESSL      -  real to complex fft :::::::::::::::
*           call drcft3(stmp,inc2r,inc3r,
*    >                  ctmp,inc2c,inc3c,
*    >                  nfft,nfft,nfft,ifwd,1.0d0,aux,naux)
*::::::::::::::::  fftpack3d -  real to complex fft :::::::::::::::
            call rc_fft3f(ctmp,
     >                     inc2c,inc3c,inc2r,inc3r,
     >	                   nfft,nfft,nfft,aux,naux)
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            do k=1,nwave
               vtmp(k) = scal1*ctmp(indx(k))
            end do
            elocal = elocal + gcdotc(nida,nidb,c1(1,n),vtmp)
         end do
      end do
      if (ispin.eq.1) elocal = 2.0d0*elocal



*     *****  get non-local psp energy *****
      call dcopy(n2ft3d*nemax,0.0d0,0,cpsi,1)
      do i=1,ni
        ia=katm(i)

*       *****   structure factor ****
        do K3=0,M3-1
          do K2=0,M2-1
            do K1=0,M1H
            CFFT(INC3C*K3+INC2C*K2+K1+1)=EX1(K1,I)*EX2(K2,I)*EX3(K3,I)
            end do
          end do
        end do
        do k=1,nwave
          exi(k)=CFFT(indx(k))
        end do

*       *****  non-local pseudopotential ****
        do L=1,LMMAX(IA)
          if(L.EQ.1 .OR. L.GE.5) then
            do K=1,NWAVE
              VTMP(K)=VNL(K,L,IA)*EXI(K)
            end do  
          else
            do K=1,NWAVE
              VTMP(K)=VNL(K,L,IA)*dcmplx(-dimag(EXI(K)),
     >                                     dble(EXI(K)))
            end do
          endif
          do N=N1(1),N2(ISPIN)
            do K=1,NWAVE
              CTMP(K)=C1(K,N)*dconjg(VTMP(K))
            end do
            SW1=gcsum(NIDA,NIDB,CTMP)/OMEGA/VNLNRM(L,IA)
            call daxpy(NWAVE2,-SW1,VTMP,1,cpsi(1,N),1)

          end do
        end do

      end do
      enlocal = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            enlocal = enlocal - gcdotc(nida,nidb,c1(1,n),cpsi(1,n))
         end do
      end do
      if (ispin.eq.1) enlocal = 2.0d0*enlocal

      

*:::::::::::::::::::::::  total energy  :::::::::::::::::::::::::::::::
      EOLD=E(1)
      E(1)=EORBIT+EION+EXC-EHARTR-PXC
      E(2)=EORBIT
      E(3)=EHARTR
      E(4)=EXC
      E(5)=EION
      E(6)=eke
      E(7)=elocal
      E(8)=enlocal
      E(9)=2.0d0*EHARTR
      E(10)=pxc

*::::::::::::::::::::::  evaluation of convergence  :::::::::::::::::::
      DELTAE=(E(1)-EOLD)/(DT*MXTIME)
      DO 690 N=1,N2(ISPIN)
        DO 690 K=1,NWAVE
          CPSI(K,N)=C2(K,N)-C1(K,N)
  690 CONTINUE
      IF(MOVE) THEN
        DO 695 I=1,NI
          FI(I,1)=(R2(I,1)-R1(I,1))/DTI(I)
          FI(I,2)=(R2(I,2)-R1(I,2))/DTI(I)
          FI(I,3)=(R2(I,3)-R1(I,3))/DTI(I)
  695   CONTINUE
      ENDIF
      DELTAC=0.0d0
      DO 720 N=N1(1),N2(ISPIN)
        ADIFF=GCDOTC(NIDA,NIDB,CPSI(1,N),CPSI(1,N))
        IF(ADIFF.GT.DELTAC) DELTAC=ADIFF
  720 CONTINUE
      DELTAC=DELTAC/DTE
      DELTAR=0.0d0
      if (move) then
         do I=1,NI
           ADIFF=dsqrt(FI(I,1)**2+FI(I,2)**2+FI(I,3)**2)
           IF(ADIFF.GT.DELTAR) DELTAR=ADIFF
         end do
      end if

      IERR=0
      RETURN
      END
