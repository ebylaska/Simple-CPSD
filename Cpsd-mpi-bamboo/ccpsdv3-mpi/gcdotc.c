
#include	"float.h"

/*************************************************
 *                                               *
 *                  gcdotc                       *
 *                                               *
 *************************************************/

/* Computes the inner produc of complex vectors in
 reciprocal space assuming that the vectors are
 defined by a real to complex 3d fft.
*/


REAL gcdotc(const int nfft, const REAL *psi1, const REAL *psi2)
{  
   int nfft2,nfft3,nffth2;
   int one=1;

   nffth2 = (nfft+2);
   nfft2 = nfft*nfft;
   nfft3 = (nfft+2)*nfft*nfft;

   return( 2.0*edot(&nfft3,psi1,&one,psi2,&one) 
           -   edot(&nfft2,psi1,&nffth2,psi2,&nffth2) 
           -   edot(&nfft2,&psi1[1],&nffth2,&psi2[1],&nffth2) );

}
