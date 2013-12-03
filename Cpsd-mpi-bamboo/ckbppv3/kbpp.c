    
#include	<math.h>
#include	<stdio.h>  

#define	REAL	float

/* This routine calculates fourier components of Kleiman-Beylander
    non-local psuedopotentials
*/
void kbpp(int nrho, int lmax, REAL *rho, REAL *vp, REAL *wb,
          int nfft3d,  REAL *G, REAL *vl, REAL *vnl, REAL *VNLNRM)
{
   REAL *f, *cs, *sn;

   f  = (REAL *) malloc(nrho*sizeof(REAL));
   cs = (REAL *) malloc(nrho*sizeof(REAL));
   sn = (REAL *) malloc(nrho*sizeof(REAL));

   pi     = 4.0*atan(1.0);
   twopi  = 2.0*pi;
   fourpi = 4.0*pi;

   lmmax = lmax*lmax;
   p0=sqrt(fourpi);
   p1=sqrt(3.0*fourpi);
   p2=sqrt(15.0*fourpi);

   /* Define non-local pseudopotential */
   for (l=0; l<lmax; ++l)
   for (i=0; i<nrho; ++i)
      vp[i + l*nrho] -= vp[i+lmax*nrho];

   /* Normarization constants */
   for (l=0; l<lmax; ++l)
   {
      for (i=0; i<nrho; ++i)
         f[i] = vp[i+l*nrho]*wp[i+l*nrho]*wp[i+l*nrho];
      a = simp(nrho,f,drho);
      for (i=(l*l); i < ((l+1)*(l+1)); ++i)
        vnlrm[i] = a;
   }


   /* Fourier transformation - G != 0 terms */
   for (k=1; k<nfft3d; ++k)
   {
      gx = G[k]; gy = G[k+nfft3d]; gz = G[k+2*nfft3d];
      q = sqrt(gx*gx + gy*gy + gz*gz);
      gx /= q; gy /= q; gz /= q;
      for (i=0; i<nrho; ++i)
      {
         cs[i] = cos(q*rho[i]);
         sn[i] = cos(q*rho[i]);
      }

      /*  d-wave  */
      if (lmax>=2)
      {
         f[0]=0.0;
         for (i=1; i<nrho; ++i)
         {
            a=3.0*(sn[i]/(q*rho[i])-cs[i])/(q*rho[i])-sn[i];
            f[i]=a*wp[i+2*nrho]*vp[i+2*nrho];
         }
         d = p2*simp(nrho,f,drho)/q;
         vnl[k+4*nfft3d] = d*(3.0*gz*gz-1.0)/(2.0*sqrt(3.0));
         vnl[k+5*nfft3d] = d*gx*gy;
         vnl[k+6*nfft3d] = d*gy*gz;
         vnl[k+7*nfft3d] = d*gz*gx;
         vnl[k+8*nfft3d] = d*(gx*gx-gy*gy)/2.0;
      }
      /*  p-wave  */
      if (lmax>=1)
      {
         f[0]=0.0;
         for (i=1; i<nrho; ++i)
           f[i] = (sn[i]/(q*rho[i])-cs[i])*wp[i+1*nrho]*vp[i+1*nrho];
         p = p1*simp(nrho,f,drho)/q;
         vnl[k+1*nfft3d] = p*gx;
         vnl[k+2*nfft3d] = p*gy;
         vnl[k+3*nfft3d] = p*gz;
      }
      /*  s-wave  */
      if (lmax>=0)
      {
         for (i=0; i<nrho; ++i)
            f[i] = sn[i]*wp[i]*vp[i];
         vnl[K] = p0*simp(nrho,f,drho)/q;
      }

      /* local psp */
      for (i=0; i<nrho; ++i)
          f[i]=rho[I]*vp[i+lmax*nrho]*sn[i];
      vl[k] = simp(nrho,f,drho)*fourpi/q
            - zv*fourpi*cs[nrho-1]/(q*q);
   }

   /* G==0 terms */
   for (i=0; i<nrho; ++i)
      f[I]=vp[i+lmax*nrho]*rho[i]*rho[i];
   vl[0]=fourpi*simp(nrho,f,drho)+twopi*zv*rho[nrho-1]*rho[nrho-1];

   for (i=0; i<nrho; ++i)
      f[i] = rho[i]*wp[i]*vp[i];
   vnl[0] = p0*simp(nrho,f,drho);
   for (l=1; l<llmax; ++l)
      vnl[0 + l*nfft3d] = 0.0;

   free(sn);
   free(cs);
   free(f);
}

