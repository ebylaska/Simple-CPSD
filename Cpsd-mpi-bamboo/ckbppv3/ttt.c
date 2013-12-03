    
#include	<math.h>
#include	<stdio.h>  

      PROGRAM MAIN
*     A driver of core routine KBPP
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      INCLUDE 'aimd.h'


      CHARACTER*2 ATOM
      CHARACTER*20 PSP

      DIMENSION VP(NRMAX,0:4),WP(NRMAX,0:4),RHO(NRMAX)
      DIMENSION VL(NWAVE),VNL(NWAVE,9),VNLNRM(9)
      DIMENSION UNITG(3,3),UNITA(3,3),G(NWAVE,3)
      DIMENSION INDX(NWAVE)
      DIMENSION RC(0:4)

      common / HAMANN / vp,wp,rho,drho,nrho,lmax,locp
      common / PSEUDO / vl,vnl,vnlnrm,zv,lmmax 
      common / INDEX1  / indx,nida,nidb
      common / LATTIC / g,unita,unitg,unit,omega,icube



      IF (IARGC().LE.0) THEN
         WRITE(*,'(A,$)') "enter pseudopotential filename =>"
         READ(*,'(A)') PSP
      ELSE
         CALL GETARG(1,PSP)
      ENDIF
*     parameters
      OPEN(9,FILE='aimd.param')
      READ(9,*,ERR=9110,END=9111) ICUBE,UNIT,IFFT
      CLOSE(9)
      IF(ICUBE.NE.JCUBE) THEN
        IERR=1
        GO TO 9999
      ENDIF
      IF(IFFT.NE.NFFT) THEN
        IERR=2
        GO TO 9999
      ENDIF

  
*     pseudopotential data
      OPEN(UNIT=11,FILE=PSP,STATUS='OLD',FORM='FORMATTED')
      READ(11,'(A2)',ERR=9116,END=9117) ATOM
      READ(11,*,ERR=9116,END=9117) ZV,AMASS,LMAX
      READ(11,*,ERR=9116,END=9117) (RC(I),I=0,LMAX)
      READ(11,*,ERR=9116,END=9117) NRHO,DRHO
      IF(NRHO.GT.NRMAX) THEN
        IERR=3
        GO TO 9999
      ENDIF
      READ(11,*,ERR=9116,END=9117)(RHO(I),(VP(I,J),J=0,LMAX),I=1,NRHO)
      READ(11,*,ERR=9116,END=9117)(RHO(I),(WP(I,J),J=0,LMAX),I=1,NRHO)
      CLOSE(11)
      WRITE(*,'(A,I1)') 'Lmax = ',LMAX
      WRITE(*,'(a,$)') "enter the highest L desired =>"
      READ(*,*) LMAX0
      LMAX=MIN(LMAX,LMAX0)

*     preparation of constants
      CALL SETUP(IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+100
        GO TO 9999
      ENDIF

      CALL KBPP(IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+200
        GO TO 9999
      ENDIF
     
      OPEN(UNIT=12,FILE='DATAOUT',STATUS='NEW',FORM='UNFORMATTED')
      WRITE(12) JCUBE,NFFT,NWAVE,UNIT
      WRITE(12) ATOM,AMASS,ZV,LMAX
      WRITE(12) (RC(I),I=0,LMAX)
      WRITE(12) (VNLNRM(I),I=1,LMMAX)
      WRITE(12) (VL(K),K=1,NWAVE)
      DO 200 L=1,LMMAX
        WRITE(12) (VNL(K,L),K=1,NWAVE)
  200 CONTINUE
      CLOSE(12)
      IERR=0
      GO TO 9999

 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 9116 IERR=16
      GO TO 9999
 9117 IERR=17
      GO TO 9999

 9999 IF(IERR.EQ.0) THEN
        write(6,*) '  VL(G=0)     =',VL(1)
        write(6,*) '  VNL(G=0,l=0)=',VNL(1,1)
        WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      ELSE
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      ENDIF
      STOP
      END

      SUBROUTINE SETUP(IERR)
*     preparation routine of KBPP
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      INCLUDE 'aimd.h'

      DIMENSION INDX(NWAVE)
      DIMENSION UNITA(3,3),UNITG(3,3),G(NWAVE,3)

      common / INDEX1  / indx,nida,nidb
      common / LATTIC / g,unita,unitg,unit,omega,icube


      call get_cube(ICUBE,UNIT,OMEGA,UNITA,UNITG)



*     index vectors for gathering and scattering
      call get_g2(NFFT,UNITG,ECUT,
     >            NIDA,NIDB,INDX,NWAVE,G)


      IERR=0
      RETURN
      END




      SUBROUTINE KBPP(IERR)
*     This routine calculates fourier components of Kleiman-Beylander
*     non-local psuedopotentials
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      INCLUDE 'aimd.h'


      DIMENSION VP(NRMAX,0:4),WP(NRMAX,0:4),RHO(NRMAX)
      DIMENSION VL(NWAVE),VNL(NWAVE,9),VNLNRM(9)
      DIMENSION UNITA(3,3),UNITG(3,3),G(NWAVE,3)
      DIMENSION F(NRMAX),CS(NRMAX),SN(NRMAX)

 
      common / HAMANN / vp,wp,rho,drho,nrho,lmax,locp
      common / PSEUDO / vl,vnl,vnlnrm,zv,lmmax 
      common / LATTIC / g,unita,unitg,unit,omega,icube


      PI=4.0d0*DATAN(1.0d0)
      TWOPI=2.0d0*PI
      FORPI=4.0d0*PI

*     Total number of non-local pseudopotentials
      LMMAX=LMAX**2

      IF(LMMAX.GT.9) THEN
        IERR=1
        RETURN
      ENDIF
      IF(NRHO.GT.NRMAX) THEN
        IERR=3
        RETURN
      ENDIF
      IF((NRHO/2)*2.EQ.NRHO) THEN
        IERR=4
        RETURN
      ENDIF

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)

*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      DO 100 L=0,LMAX-1
        DO 100 I=1,NRHO
          VP(I,L)=VP(I,L)-VP(I,LMAX)
  100 CONTINUE

*:::::::::::::::::::::  Normarization constants  ::::::::::::::::::::::
      DO 130 L=0,LMAX-1
        DO 110 I=1,NRHO
          F(I)=VP(I,L)*WP(I,L)**2
  110   CONTINUE
        A=SIMP(NRHO,F,DRHO)
        DO 120 I=L**2+1,(L+1)**2
          VNLNRM(I)=A
  120   CONTINUE
  130 CONTINUE

*======================  Fourier transformation  ======================
      DO 700 K=2,NWAVE
        Q=SQRT(G(K,1)**2+G(K,2)**2+G(K,3)**2)
        GX=G(K,1)/Q
        GY=G(K,2)/Q
        GZ=G(K,3)/Q
        DO 200 I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
  200   CONTINUE

        GO TO (600,500,400,300), LMAX+1

*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        F(1)=0.0d0
        DO 310 I=2,NRHO
          A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
          F(I)=A*WP(I,2)*VP(I,2)
  310   CONTINUE
        D=P2*SIMP(NRHO,F,DRHO)/Q
        VNL(K,5)=D*(3.0d0*GZ*GZ-1.0d0)/(2.0d0*dsqrt(3.0d0))
        VNL(K,6)=D*GX*GY
        VNL(K,7)=D*GY*GZ
        VNL(K,8)=D*GZ*GX
        VNL(K,9)=D*(GX*GX-GY*GY)/(2.0d0)

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
         F(1)=0.0d0
         DO 410 I=2,NRHO
           F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,1)*VP(I,1)
  410    CONTINUE
         P=P1*SIMP(NRHO,F,DRHO)/Q
         VNL(K,2)=P*GX
         VNL(K,3)=P*GY
         VNL(K,4)=P*GZ

*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        DO 510 I=1,NRHO
          F(I)=SN(I)*WP(I,0)*VP(I,0)
  510   CONTINUE
        VNL(K,1)=P0*SIMP(NRHO,F,DRHO)/Q

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        DO 610 I=1,NRHO
          F(I)=RHO(I)*VP(I,LMAX)*SN(I)
  610   CONTINUE
        VL(K)=SIMP(NRHO,F,DRHO)*FORPI/Q-ZV*FORPI*CS(NRHO)/(Q*Q)

  700 CONTINUE

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      DO 800 I=1,NRHO
        F(I)=VP(I,LMAX)*RHO(I)**2
  800 CONTINUE
      VL(1)=FORPI*SIMP(NRHO,F,DRHO)+TWOPI*ZV*RHO(NRHO)**2

      DO 810 I=1,NRHO
        F(I)=RHO(I)*WP(I,0)*VP(I,0)
  810 CONTINUE
      VNL(1,1)=P0*SIMP(NRHO,F,DRHO)

      DO 820 L=2,LMMAX
        VNL(1,L)=0.0d0
  820 CONTINUE

      IERR=0
      RETURN
      END

      double precision FUNCTION SIMP(N,Y,H)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      DIMENSION Y(N)

      NE=N/2
      NO=NE+1
      S=2.0d0*DSUM(NO,Y(1),2) + 4.0d0*DSUM(NE,Y(2),2)-Y(1)-Y(N)
      SIMP=S*H/3.0d0
      RETURN
      END



      subroutine get_cube(type,unit,volume,unita,unitg)

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors both in coordination           *
*     space and in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             type --- type of cube (1=SC, 2=FCC, 3=BCC)                     *
*             unit --- lattice constant                                      *
*                                                                            *
*     Outputs:                                                               *
*             volume --- volume of primitive cell                            *
*             unita  --- primitive vectors in coordination space             *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  SSCAL from BLAS                                              *
*                                                                            *
*     Last modification:  7/03/93  by R. Kawai                               *
*                                                                            *
******************************************************************************

      implicit none

*     ------------------
*     argument variables
*     ------------------
      integer type
      double precision unit, volume
      double precision unita(3,3), unitg(3,3)

*     ---------------
*     local variables
*     ---------------
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)

*     ---------------------------------------
*     primitive vectors in coordination space
*     ---------------------------------------

      if(type .eq. 1) then

*        ------------
*        simple cubic
*        ------------
         unita(1,1) = 1.0d0
         unita(2,1) = 0.0d0
         unita(3,1) = 0.0d0
         unita(1,2) = 0.0d0
         unita(2,2) = 1.0d0
         unita(3,2) = 0.0d0
         unita(1,3) = 0.0d0
         unita(2,3) = 0.0d0
         unita(3,3) = 1.0d0

      else if(type .eq. 2) then

*        -------------------
*        face centered cubic
*        -------------------
         unita(1,1) = +0.5d0
         unita(2,1) = +0.5d0
         unita(3,1) = +0.0d0
         unita(1,2) = +0.5d0
         unita(2,2) = +0.0d0
         unita(3,2) = +0.5d0
         unita(1,3) = +0.0d0
         unita(2,3) = +0.5d0
         unita(3,3) = +0.5d0
 
      else if(type .eq. 3) then
*        -------------------
*        body centered cubic
*        -------------------
         unita(1,1) = +0.5d0
         unita(2,1) = +0.5d0
         unita(3,1) = +0.5d0
         unita(1,2) = -0.5d0
         unita(2,2) = +0.5d0
         unita(3,2) = +0.5d0
         unita(1,3) = +0.5d0
         unita(2,3) = +0.5d0
         unita(3,3) = -0.5d0

      else
*        --------------------
*        unknown lattice type
*        --------------------
         write(*,*) ' get_cube: input error -- unknown lattice type.'
         stop

      endif

      call dscal(9,unit,unita,1)

*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume=dabs(volume)

      return
      end


      subroutine get_g2(ngp,unitg,ecut,ng1,ng2,pack,nwave,g)

*******************************************************************************
*                                                                             *
*     get_g2:     version 1.00                                            *
*                                                                             *
*     This routine computes reciprocal vectors within a sphere.               *
*                                                                             *
*     Ingputs:                                                                 *
*            ngp --- number of grind points in one dimention                   *
*            unitg(3,3) --- primitive vectors in reciprocal space             *
*                                                                             *
*     Outputs:                                                                *
*            ecut --- cut-off energy defining a sphere                        *
*            g(3,*) --- reciprocal vectors (Gx, Gy, Gz)                       *
*            ng1 --- number of grid points on the K1=0 plane                  *
*            ng2 --- number of grid points in upper semisphere (K1>0)         *
*            pack(*) --- index vector for gathering and scattering            *
*                          (grid points inside the sphere)                    *
*                                                                             *
*      Last Modification:  7/03/93  by Ryoichi Kawai.                         *
*                                                                             *
*******************************************************************************

      implicit none

*     ------------------
*     argument variables
*     ------------------
      integer ngp, ng0, ng1, ng2
      integer pack(*),nwave
      double precision unitg(3,3), g(nwave,3), ecut

*     ---------------
*     local variables
*     ---------------
      integer k1, k2, k3
      integer i1, i2, i3
      double precision gcut, gg, g1, g2, g3
      integer ngph, inc2, inc3

*     ---------
*     constants
*     ---------
      ngph = ngp/2
      inc2 = ngph+1
      inc3 = ngp*inc2

*     --------------
*     reset counters
*     --------------
      ng0= 0
      ng1= 0
      ng2= 0

*     ------------------
*     find cutoff energy
*     ------------------
      g1=unitg(1,1)*ngph
      g2=unitg(2,1)*ngph
      g3=unitg(3,1)*ngph
      gcut=g1*g1+g2*g2+g3*g3
      ecut = 0.5d0 * gcut

*     ---------
*     K=(0,0,0)
*     ---------
      ng1=ng1+1
      g(ng1,1) = 0.0d0
      g(ng1,2) = 0.0d0
      g(ng1,3) = 0.0d0
      pack(ng1) = 1

*     ---------------------------
*     K=(0,0,K3)  -NPH+1<K3<NPH-1
*     ---------------------------
      do k3=1,ngph-1
         g1=k3*unitg(1,3)
         g2=k3*unitg(2,3)
         g3=k3*unitg(3,3)
         gg=g1*g1+g2*g2+g3*g3
         if(gg.lt.gcut) then
            ng1=ng1+1
            g(ng1,1) = g1
            g(ng1,2) = g2
            g(ng1,3) = g3
            pack(ng1) = k3*inc3 + 1
            ng1=ng1+1
            g(ng1,1) = -g1
            g(ng1,2) = -g2
            g(ng1,3) = -g3
            pack(ng1) = (ngp-k3)*inc3 + 1
         else
            ng0 = ng0+2
         end if
      end do

*     -------------------------------------------
*     (0,K2,K3)  -NPH+1<K2<NPH-1, -NPH+1<K3<NPH-1
*     -------------------------------------------
      do k3 = -ngph+1, ngph-1
         do k2 =1, ngph-1
            g1=k2*unitg(1,2)+k3*unitg(1,3)
            g2=k2*unitg(2,2)+k3*unitg(2,3)
            g3=k2*unitg(3,2)+k3*unitg(3,3)
            gg=g1*g1+g2*g2+g3*g3
            if(gg.lt.gcut) then
               ng1 = ng1+1
               i3 = k3
               if (i3 .lt. 0) i3 = i3 + ngp
               pack(ng1) = i3*inc3 + k2*inc2 + 1
               g(ng1,1) = g1
               g(ng1,2) = g2
               g(ng1,3) = g3
               ng1 = ng1+1
               pack(ng1) = mod(ngp-k3,ngp)*inc3 
     >                   + mod(ngp-k2,ngp)*inc2 + 1
               g(ng1,1) = -g1
               g(ng1,2) = -g2
               g(ng1,3) = -g3
            else
               ng0 = ng0+2
            end if
         end do
      end do

*     ---------------------------------------------------------
*     (K1,K2,K3)   0<K1<NPH-1, -NPH+1<K2<NPH-1, -NPH+1<K3<NPH-1
*     ---------------------------------------------------------
      do k3 = -ngph+1, ngph-1
         do k2 = -ngph+1, ngph-1
            do k1 = 1, ngph-1
               g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
               g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
               g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)
               gg=g1*g1+g2*g2+g3*g3
               if(gg.lt.gcut) then
                  ng2 = ng2+1
                  i1 = k1
                  i2 = k2
                  i3 = k3
                  if (i2 .lt. 0) i2 = i2 + ngp
                  if (i3 .lt. 0) i3 = i3 + ngp
                  pack(ng1+ng2) = i3*inc3 + i2*inc2 + i1 + 1
                  g(ng1+ng2,1) = g1
                  g(ng1+ng2,2) = g2
                  g(ng1+ng2,3) = g3
               else
                  ng0 = ng0+1
               end if
            end do
         end do
      end do

*     --------------------------------------------
*     edge of primitive cells = outside the sphere
*     --------------------------------------------
      ng0 = ng0+2*ngp*ngp-ngph

      if (ng0+ng1+ng2 .ne. (ngph+1)*ngp*ngp) then
         write(*,*) ' program error: number of grids is inconsistent.'
         write(*,*) 'ng0=',ng0
         write(*,*) 'ng1=',ng1
         write(*,*) 'ng2=',ng2
         write(*,*) 'ng1+ng2=',ng1+ng2 
         write(*,*) 'ng0+ng1+ng2=',ng0+ng1+ng2,' (',(ngph+1)*ngp*ngp,')'
         stop
      end if

      return
      end
