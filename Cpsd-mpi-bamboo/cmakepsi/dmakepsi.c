#include	<stdlib.h>
#include	<stdio.h>
//#include	<stdint.h>
#include	<math.h>

/* single precision floats */
/*
#define REALTYPE  "single precision"
#define FMT1	"%f"
#define FMT10	"%10.3f %10.3f %10.3f"
#define REAL	float
#define edot	sdot
#define escal	sscal
#define ecopy	scopy
*/


/* double precision floats */
#define REALTYPE  "double precision"
#define FMT1	"%lf"
#define FMT10	"%10.3lf %10.3lf %10.3lf"
#define REAL	double
#define edot	ddot_
#define escal	dscal_
#define ecopy	dcopy_



/* define BLAS routines */
extern REAL edot();
extern void escal();
extern void ecopy();

/* define other routines */
extern REAL gcdotc(int nfft, REAL *psi1, REAL *psi2);
extern void get_cube(int itype, REAL unit, REAL *volume, REAL *unita, REAL *unitg);


main(argc,argv)
int argc;
char *argv[];
{
   int     i,j,k,l,n,ms;
   int     icube,nfft,nfft3,ispin,ne[2],neall;
   int     k1[200],k2[200],k3[200],kk2,kk3;
   REAL  sum1,scale,unit,omega,unita[9],unitg[9];
   REAL  cc[200][2],zx[2],cx[2],sx[2];
   REAL  *psi,*psi2,*psiall;
   char *filename;
   FILE *fp;

   if (argc >= 2) 
      filename = argv[1];
   else
      filename = "ELCOUT";

   /* read in aimd.param */
   fp = fopen("aimd.param","r");
     fscanf(fp,"%d",&icube);
     fscanf(fp,FMT1,&unit);
     fscanf(fp,"%d",&nfft);
     fscanf(fp,"%d %d %d",&ispin,&(ne[0]),&(ne[1]));
   fclose(fp);

   nfft3=(nfft+2)*nfft*nfft;
   neall=ne[0]+ne[1];

   psiall = (REAL *) malloc(neall*nfft3*sizeof(REAL));

   /* preparation of index vector */
   get_cube(icube,unit,&omega,unita,unitg);

   printf("\ngenerating initial psi vectors\n");
   printf(REALTYPE);printf("\n");
   printf("icube = %d\n",icube);
   printf("unit  = ");printf(FMT1,unit);printf("\n");
   printf("omega = %f\n",omega);
   printf("a1 = <"); printf(FMT10,unita[0],unita[1],unita[2]); printf(">\n");
   printf("a2 = <"); printf(FMT10,unita[3],unita[4],unita[5]); printf(">\n");
   printf("a3 = <"); printf(FMT10,unita[6],unita[7],unita[8]); printf(">\n");
   printf("g1 = <"); printf(FMT10,unitg[0],unitg[1],unitg[2]); printf(">\n");
   printf("g2 = <"); printf(FMT10,unitg[3],unitg[4],unitg[5]); printf(">\n");
   printf("g3 = <"); printf(FMT10,unitg[6],unitg[7],unitg[8]); printf(">\n");
   printf("\n");

   /* read PSI0.DAT */
   scale=1.0/sqrt(2.0);
   zx[0]=1.000; zx[1]=0.000;
   sx[0]=0.000; sx[1]=scale;
   cx[0]=scale; cx[1]=0.000;

   fp = fopen("PSI0.DAT","r");
   for (n=0; n<ne[0]; ++n)
   {
      fscanf(fp,"%d %d %d %d",&i,&j,&k,&l);
      k1[n] = i;
      k2[n] = (j+nfft)%nfft;
      k3[n] = (k+nfft)%nfft;
      if (l<0)  {cc[n][0]=sx[0];cc[n][1]=sx[1];}
      if (l==0) {cc[n][0]=zx[0];cc[n][1]=zx[1];}
      if (l>0)  {cc[n][0]=cx[0];cc[n][1]=cx[1];}
   }
   fclose(fp);

   fp = fopen(filename,"wb");
   fwrite(&icube,sizeof(int),1,fp);
   fwrite(&nfft,sizeof(int),1,fp);
   fwrite(&unit,sizeof(REAL),1,fp);
   fwrite(&ispin,sizeof(int),1,fp);
   fwrite(ne,sizeof(int),2,fp);

   i = 0; j=1; scale = 0.0;
   psi = psiall;
   for (ms=0; ms<ispin; ++ms)
   for (n=0; n<ne[ms]; ++n)
   {
     ecopy(&nfft3,&scale,&i,psi,&j);
     k = 2*(k1[n] + k2[n]*(nfft/2+1) + k3[n]*(nfft/2+1)*nfft);
     psi[k]   = cc[n][0];
     psi[k+1] = cc[n][1];
     if (k1[n]==0) 
     {
        kk2 = (nfft-k2[n])%nfft;
        kk3 = (nfft-k3[n])%nfft;
        k = 2*(k1[n] + kk2*(nfft/2+1) + kk3*(nfft/2+1)*nfft);
        psi[k]   =  cc[n][0];
        psi[k+1] = -cc[n][1];
        
     }

     psi = &psi[nfft3];
   }
   fwrite(psiall,sizeof(REAL),neall*nfft3,fp);
   fclose(fp);

   printf("checking orthonormality\n");
   psi = psiall;
   for (ms=0; ms<ispin; ++ms)
   for (n=0; n<ne[ms]; ++n)
   {
      psi2 = psi;
      for (i=n; i<ne[ms]; ++i)
      {
         printf("ms=%d  n1=%d n2=%d  <n1,ms|ms,n2>=",ms,n,i);
         printf(FMT1,gcdotc(nfft,psi,psi2));
         printf("\n");
         psi2 = &psi2[nfft3];
      }
      psi = &psi[nfft3];
   }

   free(psiall);
}



/*************************************************
 *                                               *
 *                  gcdotc                       *
 *                                               *
 *************************************************/

/* Computes the inner produc of complex vectors in 
 reciprocal space assuming that the vectors are
 defined by a real to complex 3d fft.
*/


REAL gcdotc(int nfft, REAL *psi1, REAL *psi2)
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

