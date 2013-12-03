      program aimd_size

*******************************************************************************
*     aimd_size:  version 4.00  (a new version of the former cpsize)          *
*                                                                             *
*     This program computes the number of grid points in the reciporpcal      *
*     space and other parameters that ditermines the array sizes in the       *
*     AIMD simulation package.                                                *
*                                                                             *
*     Ingputs:                                                                 *
*             cube_type --- type of cubic lattice (1=SC, 2=FCC, 3=BCC)        *
*             unit -------- lattice constant (atomic unit)                    *
*             ngp --------- number of grid points in one dimension            *
*             spin -------- spin restricted (1) or unrestricted (2)           *
*             n_elc(2) ---- number of electrons                               *
*             n_elm ------- number of elements                                *
*             n_ion ------- number of ion cores                               *
*                                                                             *
*                                                                             *
*     Outputs:                                                                *
*            aimd.h ------- header file for AIMD                              *
*            aimd.param --- parameter file for AIMD                           *
*                                                                             *
*                                                                             *
*     Last Modification  7/03/93  by R. Kawai                                 *
*                                                                             *
*******************************************************************************

      implicit none

      integer cube_type, spin, n_elc(2), n_elm, n_ion
      integer n, n1, n2, k1, k2, k3
      integer ngp, ngp0, ngp1, ngp2, ngp3, ngp5, ngph
      real unit, twopi, omega, gmax, gcut, g1, g2, g3, gg
      real unita(3,3),unitg(3,3)
      logical out_of_range

      data unita, unitg / 9*0.0, 9*0.0 /

*     ----------------
*     get lattice type
*     ----------------
      out_of_range = .true.
      do while ( out_of_range )
        write(*,'(a,$)') 'lattice type...(1) SC; (2) FCC; (3) BCC => '
        read(*,*) cube_type
        out_of_range = (cube_type.lt.1) .or. (cube_type.gt.3) 
        if (out_of_range) 
     >     write(*,*) "Lattice=",cube_type," makes no sense."
      end do

*     --------------------
*     get lattice constant
*     --------------------
      out_of_range = .true.
      do while ( out_of_range )
        write(*,'(a,$)') 'enter lattice constant => '
        read(*,*) unit
        out_of_range = unit .le. 0.
        if (out_of_range) write(*,*) " It is supposed to be positive."
      end do

*     --------------
*     get # of grids
*     --------------
      out_of_range = .true.
      do while ( out_of_range )
        write(*,'(a,$)') 'number of grids in a direction => '
        read(*,*) ngp
        n = ngp
        ngp2 = 0
        do while ( mod(n,2) .eq. 0 )
          ngp2 = ngp2 + 1
          n = n / 2
        end do
        ngp3 = 0
        do while ( mod(n,3) .eq. 0 )
          ngp3 = ngp3 + 1
          n = n / 3
        end do
        ngp5 = 0
        do while ( mod(n,5) .eq. 0 )
          ngp5 = ngp5 + 1
          n = n / 5
        end do
        out_of_range = n .ne. 1
        if ( out_of_range ) 
     >    write(*,*) " It has to be 2**i * 3**j * 5**k."
      end do

*     ----------------
*     get spin options
*     ----------------
      out_of_range = .true.
      do while ( out_of_range )
      write(*,'(a,$)') 
     > 'spin mode: (1) restricted or (2) unrestricted => '
        read(*,*) spin
        out_of_range = (spin.lt.1) .or. (spin.gt.2)
        if (out_of_range) write(*,*) "Spin=",spin," makes no sense."
      end do
       
*     ------------------
*     get # of electrons
*     ------------------
      if(spin .eq. 1) then
        out_of_range = .true.
        do while ( out_of_range )
          write(*,'(a,$)') 'number of electrons (must be even) =>'
          read(*,*) n_elc(1)
          out_of_range = (n_elc(1).le.0) .or. (mod(n_elc(1),2).eq.1)
          if (out_of_range) write(*,*) "N=",n_elc(1)," makes no sense."
        end do
        n_elc(1) = n_elc(1) / 2
        n_elc(2) = 0
      else
        out_of_range = .true.
        do while ( out_of_range )
          write(*,'(a,$)') 'number of electrons (spin majority) => '
          read(*,*) n1
          write(*,'(a,$)') 'number of electrons (spin minority) => '
          read(*,*) n2
          n_elc(1) = max(n1,n2)
          n_elc(2) = min(n1,n2)
          out_of_range = (n_elc(1).le.0) .or. (n_elc(2).lt.0)
          if (out_of_range) 
     >      write(*,*) " Negative number makes no sense."
        end do
      end if 

*     -----------------
*     get # of elements
*     -----------------
      out_of_range = .true.
      do while ( out_of_range )
        write(*,'(a,$)') "number of elements =>"  
        read(*,*) n_elm
        out_of_range = n_elm.lt.0
        if (out_of_range) write(*,*) " Negative number makes no sense."
      end do
    
*     --------------
*     get # of atoms
*     --------------
      out_of_range = .true.
      do while ( out_of_range )
        write(*,'(a,$)') "number of atoms =>"  
        read(*,*) n_ion
        out_of_range = n_ion.lt.n_elm
        if (out_of_range) write(*,*) "N=",n_ion," makes no sense."
      end do

      twopi = 8. * atan(1.0)
      ngph = ngp/2

*     -----------------
*     lattice strucutre
*     -----------------
*     simple cubic
      if(cube_type.eq.1) then
        unita(1,1)=1.0
        unita(2,2)=1.0
        unita(3,3)=1.0
        write(*,'(//a,f9.2,a)') 
     >    'lattice:  simple cubic (unit=',unit,')'

*     face centered cubic
      else if(cube_type.eq.2) then
        unita(1,1)=+0.5
        unita(2,1)=+0.5
        unita(1,2)=+0.5
        unita(3,2)=+0.5
        unita(2,3)=+0.5
        unita(3,3)=+0.5
        write(*,'(//a,f9.2,a)')
     >   'lattice:  face centered cubic (unit=',unit,')'

*     body centered cubic
      else if(cube_type.eq.3) then
        unita(1,1)=+0.5
        unita(2,1)=+0.5
        unita(3,1)=+0.5
        unita(1,2)=-0.5
        unita(2,2)=+0.5
        unita(3,2)=+0.5
        unita(1,3)=+0.5
        unita(2,3)=+0.5
        unita(3,3)=-0.5
        write(*,'(//a,f9.2,a)') 
     >   'lattice:  body centered cubic (unit=',unit,')'

      endif

      call sscal(9,unit,unita,1)

*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1)=unita(2,2)*unita(3,3)-unita(3,2)*unita(2,3)
      unitg(2,1)=unita(3,2)*unita(1,3)-unita(1,2)*unita(3,3)
      unitg(3,1)=unita(1,2)*unita(2,3)-unita(2,2)*unita(1,3)
      unitg(1,2)=unita(2,3)*unita(3,1)-unita(3,3)*unita(2,1)
      unitg(2,2)=unita(3,3)*unita(1,1)-unita(1,3)*unita(3,1)
      unitg(3,2)=unita(1,3)*unita(2,1)-unita(2,3)*unita(1,1)
      unitg(1,3)=unita(2,1)*unita(3,2)-unita(3,1)*unita(2,2)
      unitg(2,3)=unita(3,1)*unita(1,2)-unita(1,1)*unita(3,2)
      unitg(3,3)=unita(1,1)*unita(2,2)-unita(2,1)*unita(1,2)
      omega=unita(1,1)*unitg(1,1)
     &     +unita(2,1)*unitg(2,1)
     &     +unita(3,1)*unitg(3,1)
      call sscal(9,twopi/omega,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      omega=abs(omega)

      write(*,*) 'volume of primitive cell =',omega

      ngp0=0
      ngp1=0
      ngp2=0
      ngp3=0
      gmax=0.

      g1=unitg(1,1)*ngph
      g2=unitg(2,1)*ngph
      g3=unitg(3,1)*ngph
      gcut=g1*g1+g2*g2+g3*g3
      print *, gcut, ngph, unitg(1,1), unitg(2,1), unitg(3,1)

*     K=(0,0,0)
      ngp1=ngp1+1
      ngp3=ngp3+1

*     K=(0,0,K3) K3<>0
      do k3=1,ngph
        if (k3.eq.ngph) then
          ngp0=ngp0+1
          ngp3=ngp3+1
        else
          ngp3=ngp3+2
          g1=k3*unitg(1,3)
          g2=k3*unitg(2,3)
          g3=k3*unitg(3,3)
          gg=g1*g1+g2*g2+g3*g3
          if(gg.lt.gcut) then
            ngp1=ngp1+2
            if(gmax.lt.gg) gmax=gg
          else
            ngp0=ngp0+2
          endif
        end if
      end do
     
*     (0,K2,K3) 
      do k3=-ngph+1,ngph
        do k2=1,ngph-1
          if (k3.eq.ngph) then
            ngp0=ngp0+2
            ngp3=ngp3+2
          else
            ngp3=ngp3+2
            g1=k2*unitg(1,2)+k3*unitg(1,3)
            g2=k2*unitg(2,2)+k3*unitg(2,3)
            g3=k2*unitg(3,2)+k3*unitg(3,3)
            gg=g1*g1+g2*g2+g3*g3
            if(gg.lt.gcut) then
              ngp1=ngp1+2
              if(gmax.lt.gg) gmax=gg
            else
              ngp0=ngp0+2
            end if
          endif
        end do
      end do
      ngp0=ngp0+ngp
      ngp3=ngp3+ngp

*     (K1,K2,K3)
      do k3=-ngph+1,ngph
        do k2=-ngph+1,ngph
          do k1=1,ngph-1
            if (k2.eq.ngph .or. k3.eq.ngph) then
              ngp0=ngp0+1
              ngp3=ngp3+1
            else
              ngp3=ngp3+1
              g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
              g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
              g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)
              gg=g1*g1+g2*g2+g3*g3
              if(gg.lt.gcut) then
                ngp2=ngp2+1
                if(gmax.lt.gg) gmax=gg
              else
                ngp0=ngp0+1
              end if
            end if
          end do
        end do
      end do
      ngp0=ngp0+ngp*ngp
      ngp3=ngp3+ngp*ngp

      write(*,'(3(a,i2))') 'grids:',ngp,'x',ngp,'x',ngp
      write(*,*) 'semi-grids:',(ngp/2+1)*ngp*ngp,ngp3
      if ((ngp/2+1)*ngp*ngp.ne.ngp3) write(*,*) ' not consistent.'
      write(*,*) 'plane wave:',ngp1+2*ngp2
      if (ngp1+ngp2+ngp0 .ne. ngp3) write(*,*) ' not consistent.'
      write(*,*) 'array size:',ngp1+ngp2,' (',ngp1,ngp2,')'
      write(*,'(a,f7.2,a)') 'cutoff energy=',gcut,'(ryd)'
      write(*,*)

      if(spin.eq.1) then
        write(*,'(2(a,i3))') 'LDA:  up=',n_elc(1),'     down=',n_elc(1)
      else
        write(*,'(2(a,i3))') 'LSD:  up=',n_elc(1),'     down=',n_elc(2)
      end if
      write(*,'(2(a,i3))') '# of elements=',n_elm,'  # of atoms=',n_ion

      open(10,file='aimd.param')
      write(10,*) cube_type,unit,ngp
      write(10,*) spin,n_elc(1),n_elc(2)
      write(10,*) n_elm, n_ion

      open(11,file='aimd.h')
      write(11,1000) cube_type,ngp,ngp1+ngp2
      write(11,1010) n_elc(1),n_elc(2)
      write(11,1020) n_ion,n_elm
      write(11,1030)
      write(11,1040)

 1000 format(7x,'PARAMETER ( JCUBE=',i1,', NFFT=',i3,', NWAVE=',i6,' )')
 1010 format(7x,'PARAMETER ( NELC1=',i3,', NELC2=',i3,' )')
 1020 format(7x,'PARAMETER ( NATMX=',i3,', NKATMX=',i1,' )')
 1030 format(7x,'PARAMETER ( NRMAX=501 )')
 1040 format(7x,'PARAMETER ( NSHELL=1 )')


      stop
      end
