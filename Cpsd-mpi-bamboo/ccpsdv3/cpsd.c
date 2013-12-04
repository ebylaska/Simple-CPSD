#include	<stdlib.h>
#include	<math.h>
#include	<stdio.h>
#include	"float.h"

#include	"dirac_exchange.h"
#include	"vosko.h"


/* expansion coefficient of the error function */
#define CERFC	1.128379167
#define B1	0.0705230784
#define B2	0.0422820123
#define B3	0.0092705272
#define B4	0.0001520143
#define B5	0.0002765672
#define B6	0.0000430638

#define ITERLMD		30
#define CONVGLMD	1e-15


/* define other routines */
extern REAL gssum(const int nfft, const REAL *xdn);
extern REAL gsdot(const int nfft, const REAL *psi1, const REAL *psi2);
extern REAL gcdotc(const int nfft, const REAL *psi1, const REAL *psi2);

extern void phafac(const int nfft,const REAL pi, const REAL *unitg,const int nion,const REAL *r1,REAL *ex1,REAL *ex2,REAL *ex3);
extern void strfac(const int nfft, const REAL *ex1, const REAL *ex2, const REAL *ex3, REAL *exi);



void cpsd(const int inner, 
          const int move,
          const int nion,
          const int ispin,
          const int *ne,
          const int nfft,
          const int nfft3d,
          const int n2ft3d,
          const REAL *unitg,
          const REAL omega,
          const REAL dt,
          const REAL fmass,
          const int *katm,
          const REAL *amass,
          const REAL *zv,
          const int  nsh,
          const REAL *rcell, 
          const REAL cewald,
          const REAL rcut,
          const int  *masker,
          const REAL *G,
          const REAL *tg,
          const REAL *vc,
          const REAL *vg,
          int *lmmax,
          REAL *vnlnrm,
          REAL *vl,
          REAL *vnl,
          REAL *r1,
          REAL *r2,
          REAL *fion,
          REAL *c1,
          REAL *c2,
          REAL *cpsi,
          REAL *dn,
          REAL *hml,
          REAL *lmd,
          REAL *cpsdwork,
          REAL *E,REAL *deltae, REAL *deltac, REAL *deltar)
{
   int i,j,k,l,ia,ii,jj,kk,ms,itime,kc,kr,ir,jr,done;
   int nion3,nn,inc2c,inc3c,inc2r,inc3r,naux;
   int one,zero,n2ft3d2,nfft2,nffth;
   int n1[2],n2[2];
   REAL *dti,*aux,*dng,*vall,*s,*vtmp,*xcp,*xce;
   REAL *ex1,*ex2,*ex3,*rex,*rex1,*qex,*exi;
   REAL *s11,*s12,*s21,*s22,*sa0,*sa1,*st1;
   REAL *pex1,*pex2,*pex3;
   REAL dx,dy,dz,x,y,z,zz,r,aa,w,erfc,rzero,rone,rmone,adiff,f;
   REAL eorbit,eion,exc,pxc,ehartr,eke,elocal,enlocal,eold;
   REAL pi,scal1,scal2,dv,dte;
   REAL sw1,sw2,sw3,cw1[2],cw2[2],cw3[2];
   REAL rtmp;
   char tch1,tch2;


   /*  CONSTANTS */
   tch1 = 'N';
   n2ft3d2 = ispin*n2ft3d;
   nfft2   = nfft*nfft;
   nffth   = nfft/2;
   nn      = n2ft3d*(ne[0]+ne[1]);
   nion3   = 3*nion;
   inc2c   = nfft/2+1;
   inc3c   = (nfft/2+1)*nfft;
   inc2r   = (nfft+2);
   inc3r   = (nfft+2)*nfft;
   one     = 1;
   zero    = 0;
   rzero   = 0.0;
   rone    = 1.0;
   rmone   = -1.0;

   pi=4.0*atan(1.0);
   n1[0] = 0;
   n2[0] = ne[0]-1;
   n1[1] = ne[0];
   n2[1] = ne[0]+ne[1]-1;
   scal1 = 1.0/((double) (nfft*nfft*nfft));
   scal2 = 1.0/omega;
   dv    = omega*scal1;
   dte   = dt/sqrt(fmass);
   jr = 0;
   if (move)
   {
      //dti = (REAL *) malloc(nion*sizeof(REAL));
      dti = &cpsdwork[jr]; jr += nion;
      for (ii=0; ii<nion; ++ii)
         dti[ii]=dt/sqrt(amass[katm[ii]]);
   }
   naux = 2*(nfft*nfft + 10) + 40*nfft;
   aux = (REAL *) malloc(naux*sizeof(REAL));

   //dng  = (REAL *) malloc(n2ft3d*sizeof(REAL));
   //vall = (REAL *) malloc(n2ft3d*sizeof(REAL));
   //vtmp = (REAL *) malloc(n2ft3d*sizeof(REAL));
   //s    = (REAL *) malloc(n2ft3d*sizeof(REAL));
   //exi  = (REAL *) malloc(n2ft3d*sizeof(REAL));
   //xcp  = (REAL *) malloc(ispin*n2ft3d*sizeof(REAL));
   //xce  = (REAL *) malloc(n2ft3d*sizeof(REAL));

   //ex1 = (REAL *) malloc(nion*(nfft+2)*sizeof(REAL));
   //ex2 = (REAL *) malloc(nion*2*nfft*sizeof(REAL));
   //ex3 = (REAL *) malloc(nion*2*nfft*sizeof(REAL));

   //s11 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //s12 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //s21 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //s22 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //st1 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //sa1 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));
   //sa0 = (REAL *) malloc(ne[0]*ne[0]*sizeof(REAL));

   //aux  = &cpsdwork[jr]; jr += naux;
   dng  = &cpsdwork[jr]; jr += n2ft3d;
   vall = &cpsdwork[jr]; jr += n2ft3d;
   vtmp = &cpsdwork[jr]; jr += n2ft3d;
   s    = &cpsdwork[jr]; jr += n2ft3d;
   exi  = &cpsdwork[jr]; jr += n2ft3d;
   xcp  = &cpsdwork[jr]; jr += ispin*n2ft3d;
   xce  = &cpsdwork[jr]; jr += n2ft3d;

   ex1  = &cpsdwork[jr]; jr += nion*(nfft+2);
   ex2  = &cpsdwork[jr]; jr += nion*2*nfft;
   ex3  = &cpsdwork[jr]; jr += nion*2*nfft;

   s11  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   s12  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   s21  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   s22  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   st1  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   sa1  = &cpsdwork[jr]; jr += ne[0]*ne[0];
   sa0  = &cpsdwork[jr]; jr += ne[0]*ne[0];

   /*::::::::::::::::::::   |  start main loop  |   :::::::::::::::::::::::: */
   itime = 0;
   while (itime < inner)
   {
      /* increment time */
      ++itime;
      ecopy(&nn,c2,&one,c1,&one);
      ecopy(&nion3,r2,&one,r1,&one);

      /* wavefunctions in the coordination spaceJ - fftpack3d -  complex to real fft */
      ecopy(&nn,c1,&one,cpsi,&one);
      for (i=0; i<(ne[0]+ne[1]); ++i)
         cr_fft3b_(&cpsi[i*n2ft3d], 
                  &inc2c,&inc3c,&inc2r,&inc3r, 
                  &nfft,&nfft,&nfft,aux,&naux);

      /* wavefunctions in the coordination spaceJ - fftpack3d -  complex to real fft */
      /* electron spin density */
      ecopy(&n2ft3d2,&rzero,&zero,dn,&one);
      for (ms=0; ms<ispin; ++ms)
      {
         for (i=n1[ms]; i<=n2[ms]; ++i)
            for (k=0; k<n2ft3d; ++k) 
               dn[k+ms*n2ft3d] += scal2*(cpsi[k+i*n2ft3d]*cpsi[k+i*n2ft3d]);

         ecopy(&nfft2,&rzero,&zero,&dn[nfft  +ms*n2ft3d],&inc2r);
         ecopy(&nfft2,&rzero,&zero,&dn[nfft+1+ms*n2ft3d],&inc2r);
      }

      /* wavefunctions in the coordination spaceJ - fftpack3d -  complex to real fft */
      /* fourier transform of the electron density - fftpack3d - real to complex fft */
      for (k=0; k<n2ft3d; ++k)
         dng[k] = (dn[k] + dn[k+(ispin-1)*n2ft3d])*scal1;

      rc_fft3f_(dng,
               &inc2c,&inc3c,&inc2r,&inc3r,
               &nfft,&nfft,&nfft,aux,&naux);

      for (k=0; k<n2ft3d; ++k) dng[k] *= masker[k];

      /* phase factor of ion positions */
      phafac(nfft,pi,unitg,nion,r1,ex1,ex2,ex3);


      /*  kinetic energy */
      for (i=0; i<(ne[0]+ne[1]); ++i)
      {
         kk=0;
         for (k=0; k<nfft3d; ++k)
         {
            c2[kk     + i*n2ft3d] = -tg[k]*c1[kk   + i*n2ft3d];
            c2[kk + 1 + i*n2ft3d] = -tg[k]*c1[kk+1 + i*n2ft3d];
            kk += 2;
         }
      }

      /* start potential energy part */
      ecopy(&n2ft3d,&rzero,&zero,s,&one);
      ecopy(&n2ft3d,&rzero,&zero,vall,&one);

      for (ii=0; ii<nion; ++ii)
      {
         ia = katm[ii];

         /* structure factor and local pseudopotential */
         strfac(nfft,&ex1[ii*(nfft+2)],&ex2[ii*2*nfft],&ex3[ii*2*nfft],exi);
         for (k=0; k<n2ft3d; ++k) exi[k] *= masker[k];
         eaxpy(&n2ft3d,&zv[ia],exi,&one,s,&one);

         rex = &vl[ia*nfft3d];
         kk  = 0;
         for (k=0; k<nfft3d; ++k)
         {
            vtmp[kk]   = exi[kk]  *rex[k];
            vtmp[kk+1] = exi[kk+1]*rex[k];
            vall[kk]   += vtmp[kk];
            vall[kk+1] += vtmp[kk+1];
            kk += 2;
         }
         if (move)
         {
            kk = 0;
            for (k=0; k<nfft3d; ++k)
            {
               xce[k] = dng[kk+1]*vtmp[kk] - dng[kk]*vtmp[kk+1];
               kk += 2;
            }
            fion[3*ii]   = gsdot(nfft,  G,          xce);
            fion[3*ii+1] = gsdot(nfft, &G[nfft3d],  xce);
            fion[3*ii+2] = gsdot(nfft, &G[2*nfft3d],xce);
         }

         /* non-local pseudopotential */
         for (l=0; l<lmmax[ia]; ++l)
         {
            rex = &vnl[ia*9*nfft3d + l*nfft3d];

            /* s or d function */
            if ((l==0) || (l>=4))
            {
               kk = 0;
               for (k=0; k<nfft3d; ++k)
               {
                  vtmp[kk]   = exi[kk]  *rex[k];
                  vtmp[kk+1] = exi[kk+1]*rex[k];
                  kk +=2;
               }
            }

            /* p  function */
            else
            {
               kk = 0;
               for (k=0; k<nfft3d; ++k)
               {
                  vtmp[kk]   = -exi[kk+1]*rex[k];
                  vtmp[kk+1] =  exi[kk]  *rex[k];
                  kk +=2;
               }
            }
            for (i=0; i<(ne[0]+ne[1]); ++i)
            {
               sw1 = -gcdotc(nfft,vtmp,&c1[i*n2ft3d])*scal2/vnlnrm[ia*9+l];
               eaxpy(&n2ft3d,&sw1,vtmp,&one,&c2[i*n2ft3d],&one);
               if (move)
               {
                  rex1 = &c1[i*n2ft3d];
                  kk = 0;
                  for (k=0; k<nfft3d; ++k)
                  {
                     xce[k] = rex1[kk]*vtmp[kk] + rex1[kk+1]*vtmp[kk+1];
                     kk += 2;
                  }
                  if (ispin==1) sw1 *= 2;
                  fion[3*ii]   -= 2*sw1*gsdot(nfft,  G,          xce);
                  fion[3*ii+1] -= 2*sw1*gsdot(nfft, &G[nfft3d],  xce);
                  fion[3*ii+2] -= 2*sw1*gsdot(nfft, &G[2*nfft3d],xce);
               }/*move*/
            }/*i*/
         }/*l*/
      }/*ii*/

     /* Ewald summation */
     if (move)
     {
        for (ii=0; ii<nion; ++ii)
        {
           ia = katm[ii];

           /* structure factor */
           strfac(nfft,&ex1[ii*(nfft+2)],&ex2[ii*2*nfft],&ex3[ii*2*nfft],exi);
           for (k=0; k<n2ft3d; ++k) exi[k] *= masker[k];

           /* kspace ewald */
           kk = 0;
           for (k=0; k<nfft3d; ++k)
           {
              vtmp[k] = (exi[kk]*s[kk+1] - exi[kk+1]*s[kk])*vg[k];
              kk += 2;
           }
           fion[3*ii]   += gsdot(nfft,G,vtmp)           *zv[ia]*scal2;
           fion[3*ii+1] += gsdot(nfft,&G[nfft3d],vtmp)  *zv[ia]*scal2;
           fion[3*ii+2] += gsdot(nfft,&G[2*nfft3d],vtmp)*zv[ia]*scal2;

           /* realspace ewald */
           for (jj=ii+1; jj<nion; ++jj)
           {
              dx=r1[3*ii]   - r1[3*jj];
              dy=r1[1+3*ii] - r1[1+3*jj];
              dz=r1[2+3*ii] - r1[2+3*jj];
              zz=zv[katm[ii]]*zv[katm[jj]];
              sw1=0.0; sw2=0.0; sw3=0.0;
              for (l=0; l<nsh; ++l)
              {
                 x=rcell[l]        +dx;
                 y=rcell[l+nsh]    +dy;
                 z=rcell[l+nsh+nsh]+dz;
                 r=sqrt(x*x+y*y+z*z);
                 w=r/rcut;
                 aa    = (1.0+w*(B1+w*(B2+w*(B3+w*(B4+w*(B5+w*B6))))));
                 erfc  = 1.0/(aa*aa*aa*aa * aa*aa*aa*aa);
                 adiff = zz*(erfc+CERFC*w*exp(-w*w))/(r*r*r);
                 sw1 += x*adiff;
                 sw2 += y*adiff;
                 sw3 += z*adiff;
              }
              fion[3*ii]   += sw1;
              fion[3*ii+1] += sw2;
              fion[3*ii+2] += sw3;
              fion[3*jj]   -= sw1;
              fion[3*jj+1] -= sw2;
              fion[3*jj+2] -= sw3;
           }/*jj*/

        }/*ii*/
     }/*move*/

 


      /* local pseudo-, hartree-, exchange-correlation potentials  */
      Dirac_Exchange(ispin,n2ft3d,dn,xcp,xce);
      Vosko(ispin,n2ft3d,dn,xcp,xce);

      kk = 0;
      for (k=0; k<nfft3d; ++k) 
      {
         vtmp[kk]   = vall[kk]*scal2   + dng[kk]  *vc[k];
         vtmp[kk+1] = vall[kk+1]*scal2 + dng[kk+1]*vc[k];
         kk += 2;
      }
      cr_fft3b_(vtmp,
                &inc2c,&inc3c,&inc2r,&inc3r, 
                &nfft,&nfft,&nfft,aux,&naux);

      for (ms=0; ms<ispin; ++ms)
         for (i=n1[ms]; i<=n2[ms]; ++i)
         {
            for (k=0; k<n2ft3d; ++k) 
               cpsi[k+i*n2ft3d] *= (vtmp[k] + xcp[k+ms*n2ft3d]);

            rc_fft3f_(&cpsi[i*n2ft3d],
                      &inc2c,&inc3c,&inc2r,&inc3r,
                      &nfft,&nfft,&nfft,aux,&naux);

            for (k=0; k<n2ft3d; ++k) 
               cpsi[k+i*n2ft3d] = -scal1*cpsi[k+i*n2ft3d]*masker[k] + c2[k+i*n2ft3d];
         }

      /* steepest descent equations */
      ecopy(&nn,c1,&one,c2,&one);
      eaxpy(&nn,&dte,cpsi,&one,c2,&one);
      if (move)
      {
         kk = 0;
         for (ii=0; ii<nion; ++ii)
         {
            r2[kk]   = r1[kk]   + dti[i]*fion[kk];
            r2[kk+1] = r1[kk+1] + dti[i]*fion[kk+1];
            r2[kk+2] = r1[kk+2] + dti[i]*fion[kk+2];
            kk += 3;
         }
      }


      /* Lagrangian multipliers */
      for (ms=0; ms<ispin; ++ms)
      {
         kk = ne[ms]*ne[ms];
         rex1 = &c1[ms*ne[0]*n2ft3d];
         rex  = &c2[ms*ne[0]*n2ft3d];
         for (i=0; i<ne[ms]; ++i)
         {
            s22[i+i*ne[ms]] = (1.0-gcdotc(nfft,&rex[i*n2ft3d],&rex[i*n2ft3d]))*0.5/dte;
            s21[i+i*ne[ms]] = (1.0-gcdotc(nfft,&rex[i*n2ft3d],&rex1[i*n2ft3d]))*0.5;
            s12[i+i*ne[ms]] = s21[i+i*ne[ms]];
            s11[i+i*ne[ms]] = -gcdotc(nfft,&rex1[i*n2ft3d],&rex1[i*n2ft3d])*0.5*dte;
            for (j=i+1; j<ne[ms]; ++j)
            {
               s22[i+j*ne[ms]] = -gcdotc(nfft,&rex[i*n2ft3d], &rex[j*n2ft3d]) *0.5/dte;
               s21[i+j*ne[ms]] = -gcdotc(nfft,&rex[i*n2ft3d], &rex1[j*n2ft3d])*0.5;
               s12[i+j*ne[ms]] = -gcdotc(nfft,&rex1[i*n2ft3d],&rex[j*n2ft3d]) *0.5;
               s11[i+j*ne[ms]] = -gcdotc(nfft,&rex1[i*n2ft3d],&rex1[j*n2ft3d])*0.5*dte;
               s22[j+i*ne[ms]] = s22[i+j*ne[ms]];
               s21[j+i*ne[ms]] = s21[i+j*ne[ms]];
               s12[j+i*ne[ms]] = s12[i+j*ne[ms]];
               s11[j+i*ne[ms]] = s11[i+j*ne[ms]];
            }
         }
       
         ii   = 0;
         done = 0;
         ecopy(&kk,s22,&one,sa0,&one);
         while ((!done) && ((ii++)<ITERLMD))
         {
            ecopy(&kk,s22,&one,sa1,&one);
            egemm(&tch1,&tch1,&ne[ms],&ne[ms],&ne[ms],
                  &rone,
                  s21, &ne[ms],
                  sa0, &ne[ms],
                  &rone,
                  sa1, &ne[ms]);
            egemm(&tch1,&tch1,&ne[ms],&ne[ms],&ne[ms],
                  &rone,
                  sa0, &ne[ms],
                  s12, &ne[ms],
                  &rone,
                  sa1, &ne[ms]);
            egemm(&tch1,&tch1,&ne[ms],&ne[ms],&ne[ms],
                  &rone,
                  s11, &ne[ms],
                  sa0, &ne[ms],
                  &rzero,
                  st1, &ne[ms]);
            egemm(&tch1,&tch1,&ne[ms],&ne[ms],&ne[ms],
                  &rone,
                  sa0, &ne[ms],
                  st1, &ne[ms],
                  &rone,
                  sa1, &ne[ms]);
            ecopy(&kk,sa1,&one,st1,&one);
            eaxpy(&kk,&rmone,sa0,&one,st1,&one);
            adiff = fabs(st1[ieamax(&kk,st1,&one)]);
            if (adiff<CONVGLMD) 
               done = 1;
            else
               ecopy(&kk,sa1,&one,sa0,&one);
         }
         if (!done) printf("ierr=10 adiff=%f\n",adiff);
         ecopy(&kk,sa1,&one,&lmd[ms*ne[0]*ne[0]],&one);
      }

      /* correction due to the constraint */
      for (ms=0; ms<ispin; ++ms)
         egemm(&tch1,&tch1,&n2ft3d,&ne[ms],&ne[ms],
               &dte,
               &c1[ms*ne[0]*n2ft3d], &n2ft3d,
               &lmd[ms*ne[0]*ne[0]], &ne[ms],
               &rone,
               &c2[ms*ne[0]*n2ft3d],&n2ft3d);


   } /*while*/


   /* total energy calculation */

   /* hamiltonian matrix and orbital energy */
   eorbit = 0.0;
   for (ms=0; ms<ispin; ++ms)
   {
      rex  =   &c1[ms*ne[0]*n2ft3d];
      rex1 = &cpsi[ms*ne[0]*n2ft3d];
      qex  =  &hml[ms*ne[0]*ne[0]];

      for (i=0;i<ne[ms]; ++i)
      {
         qex[i+i*ne[ms]] = -gcdotc(nfft,&rex[i*n2ft3d],&rex1[i*n2ft3d]);
         eorbit += hml[i+i*ne[ms]];
         for (j=i+1;j<ne[ms]; ++j)
         {
            qex[i+j*ne[ms]] = -gcdotc(nfft,&rex[i*n2ft3d],&rex1[j*n2ft3d]);
            qex[j+i*ne[ms]] = hml[i+j*ne[ms]];
         }
      }
   }
   if (ispin==1) eorbit *= 2;



   /* ion-ion interaction energy (Ewald summation) */
   kk = 0;
   for (k=0; k<nfft3d; ++k)
   {
      vtmp[k] = (s[kk]*s[kk] + s[kk+1]*s[kk+1])*vg[k];
      kk += 2;
   }
   eion = gssum(nfft,vtmp) * 0.5*scal2 + cewald;
   for (ii=0;    ii<nion; ++ii)
   for (jj=ii+1; jj<nion; ++jj)
   {
      dx=r1[3*ii]   - r1[3*jj];
      dy=r1[1+3*ii] - r1[1+3*jj];
      dz=r1[2+3*ii] - r1[2+3*jj];
      zz=zv[katm[ii]]*zv[katm[jj]];
      for (l=0; l<nsh; ++l)
      {
         x=rcell[l]        +dx;
         y=rcell[l+nsh]    +dy;
         z=rcell[l+nsh+nsh]+dz;
         r=sqrt(x*x+y*y+z*z);
         w=r/rcut;
         aa =(1.0+w*(B1+w*(B2+w*(B3+w*(B4+w*(B5+w*B6))))));
         erfc = 1.0/(aa*aa*aa*aa);
         eion += zz*(erfc*erfc*erfc*erfc)/r;
      }
   }

   /* Hartree energy */
   kk = 0;
   for (k=0; k<nfft3d; ++k)
   {
      vtmp[k] = (dng[kk]*dng[kk] + dng[kk+1]*dng[kk+1])*vc[k];
      kk += 2;
   }
   ehartr=gssum(nfft,vtmp)*omega*0.5;

   /* exchange-correlation energy */
   exc=edot(&n2ft3d,dn,&one,xce,&one);
   pxc=edot(&n2ft3d,dn,&one,xcp,&one);
   if (ispin==1)
   {
      exc += exc;
      pxc += pxc;
   }
   else
   {
      exc += edot(&n2ft3d,&dn[n2ft3d],&one,xce,&one);
      pxc += edot(&n2ft3d,&dn[n2ft3d],&one,&xcp[n2ft3d],&one);
   }
   exc *= dv;
   pxc *= dv;


   /* Kohn-Sham kinetic energy */
   eke = 0.0;
   for (i=0; i<(ne[0]+ne[1]); ++i)
   {
      kk = 0;
      rex = &c1[i*n2ft3d];
      for (k=0; k<nfft3d; ++k)
      {
         vtmp[kk]   = rex[kk]  *tg[k];
         vtmp[kk+1] = rex[kk+1]*tg[k];
         kk += 2;
      }
      eke += gcdotc(nfft,rex,vtmp);
   }
   if (ispin==1) eke = 2*eke;

   /* Kohn-Sham V_local energy */
   elocal = gcdotc(nfft,vall,dng);

   /* get non-local psp energy */
   enlocal = 0.0;
   for (ii=0; ii<nion; ++ii)
   {
      ia = katm[ii];
      strfac(nfft,&ex1[ii*(nfft+2)],&ex2[ii*2*nfft],&ex3[ii*2*nfft],exi);
      for (k=0; k<n2ft3d; ++k) exi[k] *= masker[k];

      for (l=0; l<lmmax[ia]; ++l)
      {
         rex = &vnl[ia*9*nfft3d + l*nfft3d];

         /* s or d function */
         if ((l==0) || (l>=4))
         {
            kk = 0;
            for (k=0; k<nfft3d; ++k)
            {
               vtmp[kk]   = exi[kk]  *rex[k];
               vtmp[kk+1] = exi[kk+1]*rex[k];
               kk +=2;
            }
         }
  
         /* p  function */
         else
         {
            kk = 0;
            for (k=0; k<nfft3d; ++k)
            {
               vtmp[kk]   = -exi[kk+1]*rex[k];
               vtmp[kk+1] =  exi[kk]  *rex[k];
               kk +=2;
            }
         }
         for (i=0; i<(ne[0]+ne[1]); ++i)
         {
            sw1 = gcdotc(nfft,vtmp,&c1[i*n2ft3d]);
            enlocal += sw1*sw1*scal2/vnlnrm[ia*9 +l];
         }
      }
   }
   if (ispin==1) enlocal = 2*enlocal;


   eold = E[0];
   E[0] = eorbit + eion + exc- ehartr -pxc;
   E[1] = eorbit;
   E[2] = ehartr;
   E[3] = exc;
   E[4] = eion;
   E[5] = eke;
   E[6] = elocal;
   E[7] = enlocal;
   E[8] = 2.0*ehartr;
   E[9] = pxc;
   
   *deltae = (E[0]-eold)/(dt*inner);
   for (k=0; k<nn; ++k)
      cpsi[k] = c2[k] - c1[k];

   *deltac = 0.0;
   for (i=0; i<(ne[0]+ne[1]); ++i)
   {
      adiff = gcdotc(nfft,&cpsi[i*n2ft3d],&cpsi[i*n2ft3d]);
      if (adiff > (*deltac))
         *deltac = adiff;
   }
   *deltac /= dte;

   *deltar = 0.0;
   if (move)
   {
      for (ii=0; ii<nion; ++ii)
      {
         adiff = sqrt( fion[3*ii]*fion[3*ii]
               + fion[3*ii+1]*fion[3*ii+1]
               + fion[3*ii+2]*fion[3*ii+2] );
         if (adiff> (*deltar)) *deltar = adiff;
      }
   }


   /* deallocate memory */
/*
   free(sa1);
   free(sa0);
   free(st1);
   free(s11);
   free(s12);
   free(s21);
   free(s22);

   free(ex3);
   free(ex2);
   free(ex1);
   free(exi);
   free(xce);
   free(xcp);
   free(s);
   free(vtmp);
   free(vall);
   free(dng);
   if (move) 
      free(dti);
*/
   free(aux);
}
