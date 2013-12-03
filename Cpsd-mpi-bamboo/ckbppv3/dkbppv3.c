    
#include	<stdlib.h>
#include	<stdio.h>  
#include        <math.h>

/* single precision floats */
/*
#define REALTYPE  "single precision"
#define FMT1    "%f"
#define FMT10   "%10.3f %10.3f %10.3f"
#define REAL    float
#define edot    sdot
#define escal   sscal
#define ecopy   scopy
*/

/* double precision floats */
#define REALTYPE  "double precision"
#define FMT1    "%lf"
#define FMT10   "%10.3lf %10.3lf %10.3lf"
#define REAL    double
#define edot    ddot
#define escal   dscal
#define ecopy   dcopy


/* define BLAS routines */
extern REAL edot();
extern void escal();
extern void ecopy();


/* define other routines */
extern void get_cube(int itype, REAL unit, REAL *volume, REAL *unita, REAL *unitg);
extern void get_g(int nfft,REAL *unitg, REAL *g);
extern void kbpp(REAL zv, int lmax, int nrho, REAL drho, REAL *rho, REAL *vp, REAL *wp,
          int nfft3d,  REAL *G, REAL *vl, REAL *vnl, REAL *vnlnrm);
extern REAL gssum(int nfft, REAL *xdn);
extern REAL simp(int n, REAL *Y, REAL h);


main(int argc,char *argv[])
{
   int  i,l,nrho,lmax0,lmax,lmmax;
   int  icube,nfft,nfft3d,n2ft3d;
   int  ispin,ne[2];
   REAL unit,zv,amass,drho,rc[9];
   REAL *vnlnrm,*vl,*vnl,*g;
   REAL omega,unita[9],unitg[9];
   REAL *rho,*vp,*wp;
   char atom[2];
   char *pspname;
   FILE *fp;

   printf("\ngenerating formated pseudopotential\n");
   printf(REALTYPE);printf("\n");
   pspname = argv[1];

  /* read in aimd.param */
   fp = fopen("aimd.param","r");
     fscanf(fp,"%d",&icube);
     fscanf(fp,FMT1,&unit);
     fscanf(fp,"%d",&nfft);
     fscanf(fp,"%d %d %d",&ispin,&(ne[0]),&(ne[1]));
   fclose(fp);
   nfft3d = (nfft/2+1)*nfft*nfft;
   n2ft3d = (nfft+2)*nfft*nfft;

   /* preparation of index vector */
   get_cube(icube,unit,&omega,unita,unitg);

   g  = (REAL *) malloc(3*nfft3d*sizeof(REAL));
   get_g(nfft,unitg,g);

   printf("icube=%d\n",icube);
   printf("unit =");printf(FMT1,unit);printf("\n");
   printf("nfft =%d nfft3d=%d\n",nfft,nfft3d);
   printf("omega = %f\n",omega);
   printf("a1 = <"); printf(FMT10,unita[0],unita[1],unita[2]); printf(">\n");
   printf("a2 = <"); printf(FMT10,unita[3],unita[4],unita[5]); printf(">\n");
   printf("a3 = <"); printf(FMT10,unita[6],unita[7],unita[8]); printf(">\n");
   printf("g1 = <"); printf(FMT10,unitg[0],unitg[1],unitg[2]); printf(">\n");
   printf("g2 = <"); printf(FMT10,unitg[3],unitg[4],unitg[5]); printf(">\n");
   printf("g3 = <"); printf(FMT10,unitg[6],unitg[7],unitg[8]); printf(">\n");
   printf("\n");


   /* pseudopotential data */
   fp = fopen(pspname,"r");
     fscanf(fp,"%s",atom);
     fscanf(fp,FMT1,&zv);
     fscanf(fp,FMT1,&amass);
     fscanf(fp,"%d",&lmax);
     for (l=0; l<=lmax; ++l)
        fscanf(fp,FMT1,&rc[l]);
     fscanf(fp,"%d",&nrho);
     fscanf(fp,FMT1,&drho);

     rho = (REAL *) malloc(nrho*sizeof(REAL));
     vp  = (REAL *) malloc((lmax+1)*nrho*sizeof(REAL));
     wp  = (REAL *) malloc((lmax+1)*nrho*sizeof(REAL));
     for (i=0; i<nrho; ++i)
     {
        fscanf(fp,FMT1,&rho[i]);
        for (l=0; l<=lmax; ++l)
           fscanf(fp,FMT1,&vp[i+l*nrho]);
     }
     for (i=0; i<nrho; ++i)
     {
        fscanf(fp,FMT1,&rho[i]);
        for (l=0; l<=lmax; ++l)
           fscanf(fp,FMT1,&wp[i+l*nrho]);
     }
   fclose(fp);
   lmmax = lmax*lmax;

   printf("atom =%s\n",atom);
   printf("lmax =%d lmmax=%d\n",lmax,lmmax);
   printf("zv   =");printf(FMT1,zv);printf("\n");
   printf("rc   =");
   for (l=0; l<=lmax; ++l){printf(FMT1,rc[l]); printf(" ");}
   printf("\n");
   printf("amass=");printf(FMT1,amass);printf("\n");
   printf("nrho =%d\n",nrho);
   printf("drho =");printf(FMT1,drho);printf("\n\n");
   printf("Lmax = %d\n",lmax);
   printf("enter the highest L desired =>");
   scanf("%d",&lmax0);
   if (lmax0<lmax) lmax = lmax0;

   lmmax = lmax*lmax;
   vnlnrm = (REAL *) malloc(lmmax*sizeof(REAL));
   vl     = (REAL *) malloc(nfft3d*sizeof(REAL));
   vnl    = (REAL *) malloc(lmmax*nfft3d*sizeof(REAL));

   kbpp(zv,lmax,nrho,drho,rho,vp,wp, nfft3d,g,vl,vnl,vnlnrm);

   fp = fopen("DATAOUT","wb");
   fwrite(&icube,sizeof(int),1,fp);
   fwrite(&nfft,sizeof(int),1,fp);
   fwrite(&unit,sizeof(REAL),1,fp);
   fwrite(atom,sizeof(char),2,fp);
   fwrite(&amass,sizeof(REAL),1,fp);
   fwrite(&zv,sizeof(REAL),1,fp);
   fwrite(&lmax,sizeof(int),1,fp);
   fwrite(rc,sizeof(REAL),lmax+1,fp);
   fwrite(vnlnrm,sizeof(REAL),lmmax,fp);
   fwrite(vl,sizeof(REAL),nfft3d,fp);
   fwrite(vnl,sizeof(REAL),lmmax*nfft3d,fp);
   fclose(fp);
   printf("vlsum=%le\n",gssum(nfft,vl));

   printf("\n vnlnrm      =");
   for (l=0; l<lmmax; ++l) {printf(FMT1,vnlnrm[l]);printf(" ");}
   printf("\n");
   printf(" vl(G=0)     =");printf(FMT1,vl[0]);  printf("\n");
   printf(" vnl(G=0,l=0)=");printf(FMT1,vnl[0]); printf("\n");
   free(g);
   free(rho);
   free(vp);
   free(wp);
   free(vnl);
   free(vl);
   free(vnlnrm);
}




/*************************************************
 *                                               *
 *                get_cube                       *
 *                                               *
 *************************************************/

/*  This routine computes primitive vectors both in coordination      
      space and in reciporocal space and the volume of primitive cell. 
                                                              
      Input:  itype --- type of cube (1=SC, 2=FCC, 3=BCC)     
               unit --- lattice constant                        

      Output: volume --- volume of primitive cell                  
              unita  --- primitive vectors in coordination space    
              unitg  --- primitive vectors in reciprocal space       

      Library:  DSCAL from BLAS                                        

      Last modification:  7/03/93  by R. Kawai                           
*/

void get_cube(int itype, REAL unit, REAL *volume, REAL *unita, REAL *unitg)
{
   REAL s,v;
   int nine=9;
   int one=1;
   REAL twopi = 8.0*atan(1.0);
 

   /* primitive vectors in coordination space */

   /*simple cubic*/
   if (itype==1)
   {
      unita[0] = 1.00; unita[1] = 0.00; unita[2] = 0.00;
      unita[3] = 0.00; unita[4] = 1.00; unita[5] = 0.00;
      unita[6] = 0.00; unita[7] = 0.00; unita[8] = 1.00;

   }
   /* face centered cubic */
   else if (itype==2)
   {
      unita[0] = +0.50; unita[1] = +0.50; unita[2] = +0.00;
      unita[3] = +0.50; unita[4] = +0.00; unita[5] = +0.50;
      unita[6] = +0.00; unita[7] = +0.50; unita[8] = +0.50;
 
   }
   /* body centered cubic */
   else if (itype==3)
   {
      unita[0] = +0.50; unita[1] = +0.50; unita[2] = +0.50;
      unita[3] = -0.50; unita[4] = +0.50; unita[5] = +0.50;
      unita[6] = +0.50; unita[7] = +0.50; unita[8] = -0.50;
   }
   /* unknown lattice itype */
   else
     printf("get_cube: input error -- unknown lattice itype.\n");

   escal(&nine,&unit,unita,&one);

   /* primitive vectors in the reciprocal space */
   unitg[0+0*3] = unita[1+1*3]*unita[2+2*3] - unita[2+1*3]*unita[1+2*3];
   unitg[1+0*3] = unita[2+1*3]*unita[0+2*3] - unita[0+1*3]*unita[2+2*3];
   unitg[2+0*3] = unita[0+1*3]*unita[1+2*3] - unita[1+1*3]*unita[0+2*3];
   unitg[0+1*3] = unita[1+2*3]*unita[2+0*3] - unita[2+2*3]*unita[1+0*3];
   unitg[1+1*3] = unita[2+2*3]*unita[0+0*3] - unita[0+2*3]*unita[2+0*3];
   unitg[2+1*3] = unita[0+2*3]*unita[1+0*3] - unita[1+2*3]*unita[0+0*3];
   unitg[0+2*3] = unita[1+0*3]*unita[2+1*3] - unita[2+0*3]*unita[1+1*3];
   unitg[1+2*3] = unita[2+0*3]*unita[0+1*3] - unita[0+0*3]*unita[2+1*3];
   unitg[2+2*3] = unita[0+0*3]*unita[1+1*3] - unita[1+0*3]*unita[0+1*3];

   v = unita[0]*unitg[0] + unita[1]*unitg[1] + unita[2]*unitg[2];
   s = twopi/v;
   escal(&nine,&s,unitg,&one);

   *volume=fabs(v);
}


/**************************************************
 *                                                *
 *                 get_g                          *
 *                                                *
 **************************************************/

void get_g(int nfft,REAL *unitg, REAL *g)
{
   int i,j,k,k1,k2,k3,nffth,nfft3d,indx;
   REAL gx,gy,gz;

   nfft3d = (nfft/2+1)*nfft*nfft;
   nffth = nfft/2;

   for (k3=(-nffth+1); k3<= nffth; ++k3)
   for (k2=(-nffth+1); k2<= nffth; ++k2)
   for (k1=0; k1<= nffth; ++k1)
   {
     gx = k1*unitg[0] + k2*unitg[3] + k3*unitg[6];
     gy = k1*unitg[1] + k2*unitg[4] + k3*unitg[7];
     gz = k1*unitg[2] + k2*unitg[5] + k3*unitg[8];
     i=k1; if (i<0) i += nfft;
     j=k2; if (j<0) j += nfft;
     k=k3; if (k<0) k += nfft;
     indx = i + j*(nffth+1) + k*(nffth+1)*nfft;
     g[indx]          = gx;
     g[indx+  nfft3d] = gy;
     g[indx+2*nfft3d] = gz;
   }
}



    

REAL esum(int n, REAL *sx, int incx)
{
   int  i;
   REAL stemp;

    stemp = 0.0;
    for (i=0; i<n; ++i)
      stemp += sx[i*incx];

   return (stemp);
}

REAL gssum(int nfft, REAL *xdn)
{
   int nfft2,nfft3,nffth2;

   nffth2 = (nfft/2+1);
   nfft2 = nfft*nfft;
   nfft3 = (nfft/2+1)*nfft*nfft;

   return( 2.0*esum(nfft3,xdn,1)
           -   esum(nfft2,xdn,nffth2));

}


REAL simp(int n, REAL *Y, REAL h)
{
   int i,ne,no;
   REAL s,s2,s4;
   //ne = n/2;
   //no = ne+1;

   s2 = 0.0;
   for (i=0; i<n; i+=2)  
      s2 += Y[i];
   s4 = 0.0;
   for (i=1; i<n; i+=2) 
      s4 += Y[i];

   s = 2.0*s2 + 4.0*s4 - Y[0] - Y[n-1];
   s *= (h/3.0);
   
  return s;
}

    

/* This routine calculates fourier components of Kleiman-Beylander
    non-local psuedopotentials
*/
void kbpp(REAL zv, int lmax, int nrho, REAL drho, REAL *rho, REAL *vp, REAL *wp,
          int nfft3d,  REAL *G, REAL *vl, REAL *vnl, REAL *vnlnrm)
{
   int i,l,k,lmmax;
   REAL pi,twopi,fourpi;
   REAL p0,p1,p2,p,d,gx,gy,gz,a,q;
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
        vnlnrm[i] = a;
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
         sn[i] = sin(q*rho[i]);
      }

      /*  d-wave  */
      if (lmax>2)
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
      if (lmax>1)
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
      if (lmax>0)
      {
         for (i=0; i<nrho; ++i)
            f[i] = sn[i]*wp[i]*vp[i];
         vnl[k] = p0*simp(nrho,f,drho)/q;
      }

      /* local psp */
      for (i=0; i<nrho; ++i)
          f[i]=rho[i]*vp[i+lmax*nrho]*sn[i];
      vl[k] = simp(nrho,f,drho)*fourpi/q
            - zv*fourpi*cs[nrho-1]/(q*q);
   }


   /* G==0 terms */
   for (i=0; i<nrho; ++i)
   {
      f[i] = vp[i + lmax*nrho] * (rho[i]*rho[i]);
   }
   a = simp(nrho,f,drho);
   vl[0]=fourpi*a + twopi*zv*rho[nrho-1]*rho[nrho-1];

   for (i=0; i<nrho; ++i)
      f[i] = rho[i]*wp[i]*vp[i];
   vnl[0] = p0*simp(nrho,f,drho);
   for (l=1; l<lmmax; ++l)
      vnl[0 + l*nfft3d] = 0.0;

   free(sn);
   free(cs);
   free(f);
}

