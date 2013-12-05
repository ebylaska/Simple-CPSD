    

#include	"float.h"


/**************************************************
 *                                                *
 *                 get_masker                     *
 *                                                *
 **************************************************/

void get_masker(const int nfft, const REAL ggcut, const REAL *unitg, int *masker1)
{
   int i,j,k,k1,k2,k3,nffth,nfft3d,indx;
   REAL gx,gy,gz,gg;

   nfft3d = (nfft/2+1)*nfft*nfft;
   nffth = nfft/2;
   for (k=0; k<nfft3d; ++k) masker1[k] = 0;

   for (k3=(-nffth+1); k3<nffth; ++k3)
   for (k2=(-nffth+1); k2<nffth; ++k2)
   for (k1=0; k1<nffth; ++k1)
   {
     gx = k1*unitg[0] + k2*unitg[3] + k3*unitg[6];
     gy = k1*unitg[1] + k2*unitg[4] + k3*unitg[7];
     gz = k1*unitg[2] + k2*unitg[5] + k3*unitg[8];
     gg = gx*gx + gy*gy + gz*gz;

     i=k1; if (i<0) i += nfft;
     j=k2; if (j<0) j += nfft;
     k=k3; if (k<0) k += nfft;
     indx = i + j*(nffth+1) + k*(nffth+1)*nfft;
     if (gg < ggcut) masker1[indx] = 1;

   }
}
