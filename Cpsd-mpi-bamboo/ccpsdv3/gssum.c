
#include	"float.h"





/*************************************************
 *                                               *
 *                  gssum                        *
 *                                               *
 *************************************************/

/* Computes the inner produc of complex vectors in
 reciprocal space assuming that the vectors are
 defined by a real to complex 3d fft.
*/

REAL esum(const int n, const REAL *sx, const int incx)
{     
   int  i;
   REAL stemp;
   stemp = 0.0;
   for (i=0; i<n; ++i)
      stemp += sx[i*incx];
   return (stemp);
}


REAL gssum(const int nfft, const REAL *xdn)
{  
   int nfft2,nfft3,nffth2;

   nffth2 = (nfft/2+1);
   nfft2 = nfft*nfft;
   nfft3 = (nfft/2+1)*nfft*nfft;

   return( 2.0*esum(nfft3,xdn,1)
           -   esum(nfft2,xdn,nffth2));

}
