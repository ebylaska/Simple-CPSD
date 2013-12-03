*     ***********************************
*     *					*
*     *	        cr_fft3b		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of     
*      a three dimensional complex to real
*      inverse fft                                
*           A(nx,ny,nz) <- FFT3^(-1)[A(kx,ky,kz)] 
*                                                 
*      Entry - 					  
*              tmp: tempory work space must be at 
*                    least the size of (complex)  
*                    (nfft*nfft + 1) + 10*nfft    
*                                                 
*       Exit - A is transformed and the imaginary 
*              part of A is set to zero           
*                                                 

      subroutine cr_fft3b(A,
     >                    inc2c,inc3c,
     >                    inc2r,inc3r,
     >                    nx,ny,nz,tmp,ntmp)

      implicit none
      
      real*4  A(*)
      integer inc2c,inc3c
      integer inc2r,inc3r
      integer nx,ny,nz
      real*4  tmp(*)
      integer  ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tz,ty,tx

      tz = 2*nz+1
      ty = 2*ny+1
      tx = nx+3

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,ky,nz) <- fft1d^(-1)[A(kx,ky,kz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call cffti(nz,tmp(tz))
      do j=1,ny
      do i=1,(nx/2+1)
         indx = (2*i-1) + (j-1)*(inc2r) 
         call ccopy(nz,A(indx),inc3c,tmp,1)
         call cfftb(nz,tmp,tmp(tz))
         call ccopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do


*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call cffti(ny,tmp(ty))
      do k=1,nz
      do i=1,(nx/2+1)
         indx = (2*i-1) + (k-1)*(inc3r) 
         call ccopy(ny,A(indx),inc2c,tmp,1)
         call cfftb(ny,tmp,tmp(ty))
         call ccopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)]  ***
*     ***   1-d complex to real transforms          ***
*     *************************************************
      call rffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2r) + (k-1)*(inc3r)
         call scopy((nx+2),A(indx),1,tmp,1)
         do i=2,nx
            tmp(i) = tmp(i+1)
         end do
         call rfftb(nx,tmp,tmp(tx))
         tmp(nx+1) = 0.0d0
         tmp(nx+2) = 0.0d0
         call scopy((nx+2),tmp,1,A(indx),1)
      end do
      end do
      return
      end




*     ***********************************
*     *					*
*     *	       rc_fft3f			*
*     *					*
*     ***********************************

      subroutine rc_fft3f(A,
     >                    inc2c,inc3c,
     >                    inc2r,inc3r,
     >                    nx,ny,nz,tmp,ntmp)


*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx,ny,nz)]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*                                                   *
*****************************************************
      implicit none
      real*4   A(*)
      integer inc2c,inc3c
      integer inc2r,inc3r
      integer nx,ny,nz
      real*4  tmp(*)
      integer ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tx,ty,tz

      tx = nx   + 3
      ty = 2*ny + 1
      tz = 2*nz + 1


*     ********************************************
*     ***     do fft along nx dimension        ***
*     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
*     ********************************************
      call rffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2r) + (k-1)*(inc3r)
         call scopy((nx+2),A(indx),1,tmp,1)
         call rfftf(nx,tmp,tmp(tx))
         do i=nx,2,-1
            tmp(i+1) = tmp(i)
         end do
         tmp(2)    = 0.0d0
         tmp(nx+2) = 0.0d0
         call scopy((nx+2),tmp,1,A(indx),1)
      end do
      end do

*     ********************************************
*     ***     do fft along ny dimension        ***
*     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
*     ********************************************
      call cffti(ny,tmp(ty))
      do k=1,nz
      do i=1,(nx/2+1)
         indx = (2*i-1) + (k-1)*(inc3r) 
         call ccopy(ny,A(indx),inc2c,tmp,1)
         call cfftf(ny,tmp,tmp(ty))
         call ccopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do


*     ********************************************
*     ***     do fft along nz dimension        ***
*     ***   A(kx,ky,kz) <- fft1d[A(kx,ky,nz)]  ***
*     ********************************************
      call cffti(nz,tmp(tz))
      do j=1,ny
      do i=1,(nx/2+1)
         indx = (2*i-1) + (j-1)*(inc2r) 
         call ccopy(nz,A(indx),inc3c,tmp,1)
         call cfftf(nz,tmp,tmp(tz))
         call ccopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *	        cc_fft3b		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of     
*      a three dimensional complex to complex     
*      inverse fft                                
*           A(nx,ny,nz) <- FFT3^(-1)[A(kx,ky,kz)] 
*                                                 
*      Entry - 					  
*              tmp: tempory work space must be at 
*                    least the size of (complex)  
*                    (nfft*nfft + 1) + 10*nfft    
*                                                 
*       Exit - A is transformed and the imaginary 
*              part of A is set to zero           
*                                                 

      subroutine cc_fft3b(A,inc2c,inc3c,
     >                    nx,ny,nz,tmp,ntmp)

      implicit none
      
      complex*8 A(*)
      integer    inc2c,inc3c
      integer    nx,ny,nz
      complex*8 tmp(*)
      integer    ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tz,ty,tx

      tz = nz+1
      ty = ny+1
      tx = nx+1


*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,ky,nz) <- fft1d^(-1)[A(kx,ky,kz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call cffti(nz,tmp(tz))
      do j=1,ny
      do i=1,nx
         indx = i + (j-1)*(inc2c) 
         call ccopy(nz,A(indx),inc3c,tmp,1)
         call cfftb(nz,tmp,tmp(tz))
         call ccopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do


*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call cffti(ny,tmp(ty))
      do k=1,nz
      do i=1,nx
         indx = i + (k-1)*(inc3c) 
         call ccopy(ny,A(indx),inc2c,tmp,1)
         call cfftb(ny,tmp,tmp(ty))
         call ccopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call cffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2c) + (k-1)*(inc3c)
         call ccopy(nx,A(indx),1,tmp,1)
         call cfftb(nx,tmp,tmp(tx))
         call ccopy(nx,tmp,1,A(indx),1)
      end do
      end do
      return
      end




*     ***********************************
*     *					*
*     *	       cc_fft3f			*
*     *					*
*     ***********************************

      subroutine cc_fft3f(A,inc2c,inc3c,
     >                    nx,ny,nz,tmp,ntmp)


*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx,ny,nz)]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*                                                   *
*****************************************************
      implicit none
      complex*8 A(*)
      integer    inc2c,inc3c
      integer    nx,ny,nz
      complex*8 tmp(*)
      integer    ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tx,ty,tz

      tx = nx + 1
      ty = ny + 1
      tz = nz + 1


*     ********************************************
*     ***     do fft along nx dimension        ***
*     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
*     ********************************************
      call cffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2c) + (k-1)*(inc3c)
         call ccopy(nx,A(indx),1,tmp,1)
         call cfftf(nx,tmp,tmp(tx))
         call ccopy(nx,tmp,1,A(indx),1)
      end do
      end do

*     ********************************************
*     ***     do fft along ny dimension        ***
*     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
*     ********************************************
      call cffti(ny,tmp(ty))
      do k=1,nz
      do i=1,nx
         indx = i + (k-1)*(inc3c) 
         call ccopy(ny,A(indx),inc2c,tmp,1)
         call cfftf(ny,tmp,tmp(ty))
         call ccopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do


*     ********************************************
*     ***     do fft along nz dimension        ***
*     ***   A(kx,ky,kz) <- fft1d[A(kx,ky,nz)]  ***
*     ********************************************
      call cffti(nz,tmp(tz))
      do j=1,ny
      do i=1,nx
         indx = i + (j-1)*(inc2c) 
         call ccopy(nz,A(indx),inc3c,tmp,1)
         call cfftf(nz,tmp,tmp(tz))
         call ccopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do

      return
      end

