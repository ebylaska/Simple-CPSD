#include	<math.h>
#include	"float.h"

/*********************************************
 *                                           *
 *                phafac                     *
 *                                           *
 *********************************************/

void phafac(const int nfft, const REAL pi, const REAL *unitg, const int nion, const REAL *r1, 
            REAL *ex1, REAL *ex2, REAL *ex3)
{
   int ii,k,nffth;
   REAL sw1,sw2,sw3,cw1[2],cw2[2],cw3[2];
   REAL *rex,*rex1,*qex;

    nffth = nfft/2;

   /* phase factor of ion positions */
   for (ii=0; ii<nion; ++ii)
   {
      sw1=unitg[0]*r1[3*ii]+unitg[1]*r1[3*ii+1]+unitg[2]*r1[3*ii+2]+pi;
      sw2=unitg[3]*r1[3*ii]+unitg[4]*r1[3*ii+1]+unitg[5]*r1[3*ii+2]+pi;
      sw3=unitg[6]*r1[3*ii]+unitg[7]*r1[3*ii+1]+unitg[8]*r1[3*ii+2]+pi;
      cw1[0] = cos(sw1); cw1[1] = -sin(sw1);
      cw2[0] = cos(sw2); cw2[1] = -sin(sw2);
      cw3[0] = cos(sw3); cw3[1] = -sin(sw3);

      ex1[0 + ii*(nfft+2)] = 1.0; ex1[1 + ii*(nfft+2)] = 0.0;
      ex2[0 + ii*2*nfft]   = 1.0; ex2[1 + ii*2*nfft]   = 0.0;
      ex3[0 + ii*2*nfft]   = 1.0; ex3[1 + ii*2*nfft]   = 0.0;

      rex  = &ex1[ii*(nfft+2)];
      rex1 = &rex[2];
      for (k=1; k<=nffth; ++k)
      {
         rex1[0] = rex[0]*cw1[0] - rex[1]*cw1[1];
         rex1[1] = rex[0]*cw1[1] + rex[1]*cw1[0];
         ++rex;  ++rex;
         ++rex1; ++rex1;
      }

      rex  = &ex2[ii*2*nfft];
      rex1 = &rex[2];
      qex  = &ex2[(ii+1)*2*nfft-2];
      for (k=1; k<=nffth; ++k)
      {
         rex1[0] = rex[0]*cw2[0] - rex[1]*cw2[1];
         rex1[1] = rex[0]*cw2[1] + rex[1]*cw2[0];
         qex[0] =  rex1[0];
         qex[1] = -rex1[1];
         ++rex;  ++rex;
         ++rex1; ++rex1;
         --qex;  --qex;
      }
      rex  = &ex3[ii*2*nfft];
      rex1 = &rex[2];
      qex  = &ex3[(ii+1)*2*nfft-2];
      for (k=1; k<=nffth; ++k)
      {
         rex1[0] = rex[0]*cw3[0] - rex[1]*cw3[1];
         rex1[1] = rex[0]*cw3[1] + rex[1]*cw3[0];
         qex[0] =  rex1[0];
         qex[1] = -rex1[1];
         ++rex;  ++rex;
         ++rex1; ++rex1;
         --qex;  --qex;
      }
      ex1[ii*(nfft+2)+nfft]   = 0.0;
      ex1[ii*(nfft+2)+nfft+1] = 0.0;
      ex2[ii*2*nfft+nfft]     = 0.0;
      ex2[ii*2*nfft+nfft+1]   = 0.0;
      ex3[ii*2*nfft+nfft]     = 0.0;
      ex3[ii*2*nfft+nfft+1]   = 0.0;
   }

}


/*********************************************
 *                                           *
 *                strfac                     *
 *                                           *
 *********************************************/

/*  
  This routine computes the structure factor from the phase factors.

  Entry: nfft - grid size
         ex1,ex2,ex3 - phase factors
  Exit:  exi - structure factor
*/
void	strfac(const int nfft,
               const REAL *ex1,
               const REAL *ex2,
               const REAL *ex3,
               REAL *exi)
{
   int i,j,k,ir,jr,kr,kk,nffth;

   nffth = nfft/2;
   kk = 0;
   kr = 0;
   for (k=0; k<nfft; ++k)
   {
      jr = 0;
      for (j=0; j<nfft; ++j)
      {
         ir = 0;
         for (i=0; i<=nffth; ++i)
         {
            exi[kk]   = ex1[ir]   * (ex2[jr]*ex3[kr]   - ex2[jr+1]*ex3[kr+1])
                      - ex1[ir+1] * (ex2[jr]*ex3[kr+1] + ex2[jr+1]*ex3[kr]);
            exi[kk+1] = ex1[ir+1] * (ex2[jr]*ex3[kr]   - ex2[jr+1]*ex3[kr+1])
                      + ex1[ir]   * (ex2[jr]*ex3[kr+1] + ex2[jr+1]*ex3[kr]);
            kk += 2;
            ir +=2;
         }
         jr += 2;
      }
      kr +=2;
   }
}

