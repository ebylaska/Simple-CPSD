#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>

#include	"float.h"
#include	"Parallel.h"
#include	"d3db.h"

/* define other routines */
extern void get_cube(int itype, REAL unit, REAL *volume, REAL *unita, REAL *unitg);
extern void get_g(int nfft,REAL *unitg, REAL *g);
extern void get_masker(const int nfft, const REAL ggcut, const REAL *unitg, int *masker1);
extern REAL gcdotc(const int nfft, const REAL *psi1, const REAL *psi2);
extern REAL gssum(const int nfft, const REAL *xdn);
extern REAL current_second();
extern void eigsrt(REAL *D, REAL *V, int n);


extern void cpsd(const int inner, 
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
          REAL *E,REAL *deltae, REAL *deltac, REAL *deltar);


main(int argc, char *argv[])
{
   int i,j,k,kk,l,n,ia,ii,ms,nn,one,zero,oprint,mapping;
   int i1,i2,i3,nsh;
   int nfft,nfft3d,n2ft3d;
   int ispin,imove,icube,inner,outer,ncut,icount,done;
   int nkatm,nion;
   int lmax0,icube0,nfft0,lwork,info;
   int *lmax,*lmmax,*natm,*katm,*masker,*masker1;
   int ne[2],n1[2],n2[2];

   REAL *r1,*r2,*fion;
   REAL *vl, *vnl, *vnlnrm, *zv,*amass;
   REAL **rc;
   REAL *eig,*c1,*c2,*cpsi,*lmd,*hml,*psi,*psi2;
   REAL *work;
   REAL *dn,*xcp,*xce;
   REAL *tg,*vc,*vg,*G,*rcell;

   int  cpsdlwork;
   REAL *cpsdwork;

   REAL unit,unit0,dt,fmass,tole,tolc,tolr,rzero,rone;
   REAL rcut,pi,fourpi,alpha,w,gg,gx,gy,gz,ecut,gcut;
   REAL cpu1,cpu2,cpu3,cpu4,cpu5,t1,t2,t3,t4,av;
   REAL q,aa,rs,zz,z,cewald;
   REAL am,tmp,omega,unita[9],unitg[9],gc[3],cm[3];
   REAL E[10],deltae,deltar,deltac,ev;

   char **fname,**atom;
   char tch1,tch2,tch3;
   FILE *fp;


/*****************************|  PROLOGUE  |****************************/


   cpu1 = current_second();
   one   = 1;
   zero  = 0;
   rone  = 1.0;
   rzero = 0.0;

   Parallel_init(&argc,&argv);
   oprint = (Parallel_taskid()==0);

   if (oprint)
   {
      printf("          ****************************************************\n");
      printf("          *                                                  *\n");
      printf("          *     Car-Parrinello microcluster calculation      *\n");
      printf("          *                                                  *\n");
      printf("          *         [ steepest descent minimization ]        *\n");
      printf("          *         [      ");
      printf(REALTYPE);
      printf("         ]        *\n");
      printf("          *                                                  *\n");
      printf("          *         [   parallel MPI implementation ]        *\n");
      printf("          *            version #3.00   12/05/2013            *\n");
      printf("          *                                                  *\n");
      printf("          ****************************************************\n");
      message(1);
   }

   /* read parameters and flags */
   if (Parallel_taskid()==0)
   {
      fp = fopen("CONTROL","r");
      fscanf(fp,"%d %d", &ispin, &imove);
      fscanf(fp,"%d", &icube);
      fscanf(fp,FMT1,&unit);
      fscanf(fp,FMT1,&dt);
      fscanf(fp,FMT1,&fmass);
      fscanf(fp,"%d %d", &inner, &outer);
      fscanf(fp,FMT1,&tole);
      fscanf(fp,FMT1,&tolc);
      fscanf(fp,FMT1,&tolr);
      if (fscanf(fp,"%d",&mapping)==EOF) 
         mapping = 1;
      fclose(fp);
   }
   Parallel_ibcast(0,1,&ispin);
   Parallel_ibcast(0,1,&imove);
   Parallel_ibcast(0,1,&icube);
   Parallel_ibcast(0,1,&inner);
   Parallel_ibcast(0,1,&outer);
   Parallel_rbcast(0,1,&unit);
   Parallel_rbcast(0,1,&dt);
   Parallel_rbcast(0,1,&tole);
   Parallel_rbcast(0,1,&tolc);
   Parallel_rbcast(0,1,&tolr);
   Parallel_ibcast(0,1,&mapping);

   /* read in ELCIN header */
   if (Parallel_taskid()==0)
   {
      fp = fopen("ELCIN","rb");
      fread(&icube,sizeof(int),1,fp);
      fread(&nfft,sizeof(int),1,fp);
      fread(&unit,sizeof(REAL),1,fp);
      fread(&ispin,sizeof(int),1,fp);
      fread(ne,sizeof(int),2,fp);
   }
   Parallel_ibcast(0,1,&icube);
   Parallel_ibcast(0,1,&nfft);
   Parallel_rbcast(0,1,&unit);
   Parallel_ibcast(0,1,&ispin);
   Parallel_ibcast(0,2,ne);

   /* initialize d3db */
   d3db_init(nfft,nfft,nfft,mapping);
   nfft3d = d3db_nfft3d();
   n2ft3d = d3db_n2ft3d();

   /* allocate electronic data */
   //nfft3d = (nfft/2+1)*nfft*nfft;
   //n2ft3d = (nfft+2)*nfft*nfft;
   c1   = (REAL *) malloc((ne[0]+ne[1])*n2ft3d*sizeof(REAL));
   c2   = (REAL *) malloc((ne[0]+ne[1])*n2ft3d*sizeof(REAL));
   cpsi = (REAL *) malloc((ne[0]+ne[1])*n2ft3d*sizeof(REAL));
   dn  = (REAL *) malloc(ispin*n2ft3d*sizeof(REAL));
   xcp = (REAL *) malloc(2*n2ft3d*sizeof(REAL));
   xce = (REAL *) malloc(2*n2ft3d*sizeof(REAL));
   hml = (REAL *) malloc((ne[0]*ne[0]+ne[1]*ne[1])*sizeof(REAL));
   lmd = (REAL *) malloc((ne[0]*ne[0]+ne[1]*ne[1])*sizeof(REAL));
   eig = (REAL *) malloc((ne[0]+ne[1])*sizeof(REAL));
   lwork = 4*ne[0];
   work = (REAL *) malloc(lwork*sizeof(REAL));

   for (i=0; i<(ne[0]+ne[1]); ++i)
      d3db_c_read(fp,&c2[i*n2ft3d],c1,cpsi);

   nn = (ne[0]+ne[1])*n2ft3d;
   //fread(c2,sizeof(REAL),nn,fp);
   ecopy(&nn,c2,&one,c1,&one);

   if (Parallel_taskid()==0)
      fclose(fp);

   /* define lattice */
   get_cube(icube,unit,&omega,unita,unitg);

   pi = 4.0*atan(1.0);
   fourpi = 4.0*pi;
   alpha = 0.0;
   if (icube==1) alpha = 1.76011888;
   if (icube==2) alpha = 1.791753;
   if (icube==3) alpha = 1.791860;
   rcut = unit/pi;
   ncut = 1;

   /* preparation of index vector */
   G  = (REAL *) malloc(3*nfft3d*sizeof(REAL));
   get_g(nfft,unitg,G);

   /* find cutoff energy */
   gx=unitg[0]*(nfft/2);
   gy=unitg[1]*(nfft/2);
   gz=unitg[2]*(nfft/2);
   gcut=gx*gx + gy*gy + gz*gz;
   ecut=0.50*gcut;


   /* set kinetic, coulomb and ewald vectors */
   tg  = (REAL *) malloc(nfft3d*sizeof(REAL));
   vc  = (REAL *) malloc(nfft3d*sizeof(REAL));
   vg  = (REAL *) malloc(nfft3d*sizeof(REAL));
   w = 0.25*rcut*rcut;
   for (k=1; k<nfft3d; ++k)
   {
      gg    = G[k]*G[k]+G[k+nfft3d]*G[k+nfft3d]+G[k+2*nfft3d]*G[k+2*nfft3d];
      tg[k] = 0.5*gg;
      vc[k] = fourpi/gg;
      vg[k] = vc[k]*exp(-w*gg);
   }
   tg[0] = 0.0; 
   vc[0] = 0.0; 
   vg[0] = 0.0;


   /* define the mask array */
   masker1= (int *) malloc(nfft3d*sizeof(int));
   masker = (int *) malloc(n2ft3d*sizeof(int));
   get_masker(nfft,gcut,unitg,masker1);
   kk = 0;
   for (k=0; k<nfft3d; ++k)
   {
      masker[kk]   = masker1[k];
      masker[kk+1] = masker1[k];
      kk += 2;
   }
   for (k=0; k<nfft3d; ++k) tg[k] *= masker1[k];
   for (k=0; k<nfft3d; ++k) vc[k] *= masker1[k];
   for (k=0; k<nfft3d; ++k) vg[k] *= masker1[k];

   /*  read ionic structure */
   if (Parallel_taskid()==0)
   {
      fp = fopen("IONIN","r");
      fscanf(fp,"%d",&nkatm);
   }
   Parallel_ibcast(0,1,&nkatm);

   lmax   = (int *)   malloc(nkatm*sizeof(int));
   lmmax  = (int *)   malloc(nkatm*sizeof(int));
   natm   = (int *)   malloc(nkatm*sizeof(int));
   fname  = (char **) malloc(nkatm*sizeof(char *));
   atom   = (char **) malloc(nkatm*sizeof(char *));
   for (ia=0; ia<nkatm; ++ia)
   {
      natm[ia]  = 0;
      fname[ia] = (char *) malloc(8*sizeof(char));
      atom[ia]  = (char *) malloc(2*sizeof(char));
   }

   if (Parallel_taskid()==0)
   {
     for (ia=0; ia<nkatm; ++ia)
     {
        fscanf(fp,"%1s%s%1s %d",&tch1,fname[ia],&tch2,&lmax[ia]);
        lmmax[ia] = lmax[ia]*lmax[ia];
     }
     fscanf(fp,"%d",&nion);
   }
   Parallel_ibcast(0,nkatm,lmax);
   Parallel_ibcast(0,nkatm,lmmax);
   Parallel_ibcast(0,1,&nion);

   katm = (int *)  malloc(nion*sizeof(int));
   r2   = (REAL *) malloc(3*nion*sizeof(REAL));
   r1   = (REAL *) malloc(3*nion*sizeof(REAL));
   fion = (REAL *) malloc(3*nion*sizeof(REAL));

   if (Parallel_taskid()==0)
   {
      for (ii=0; ii<nion; ++ii)
      {
         fscanf(fp,"%d",&katm[ii]); katm[ii] -= 1;
         fscanf(fp,FMT1,&r2[0+ii*3]); fscanf(fp,FMT1,&r2[1+ii*3]); fscanf(fp,FMT1,&r2[2+ii*3]);
         fscanf(fp,FMT1,&r1[0+ii*3]); fscanf(fp,FMT1,&r1[1+ii*3]); fscanf(fp,FMT1,&r1[2+ii*3]);
         natm[katm[ii]] += 1;
      }
      fclose(fp);
   }
   Parallel_ibcast(0,nion,katm);
   Parallel_rbcast(0,3*nion,r2);
   nn = 3*nion;
   ecopy(&nn,r2,&one,r1,&one);

   /* read pseudopotential */
   amass = (REAL *) malloc(nkatm*sizeof(REAL));
   zv    = (REAL *) malloc(nkatm*sizeof(REAL));

   rc  = (REAL **) malloc(nkatm*sizeof(REAL *));
   for (ia=0; ia<nkatm; ++ia) rc[ia] = (REAL *) malloc((lmax[ia]+1)*sizeof(REAL));

   vnlnrm  = (REAL *) malloc(nkatm*9*sizeof(REAL));
   vl      = (REAL *) malloc(nkatm*nfft3d*sizeof(REAL));
   vnl     = (REAL *) malloc(nkatm*nfft3d*9*sizeof(REAL));

   for (ia=0; ia<nkatm; ++ia)
   {
      if (Parallel_taskid()==0)
      {
         fp = fopen(fname[ia],"rb");
         fread(&icube0,sizeof(int),1,fp);
         fread(&nfft0,sizeof(int),1,fp);
         fread(&unit0,sizeof(REAL),1,fp);
         fread(atom[ia],sizeof(char),2,fp);
         fread(&amass[ia],sizeof(REAL),1,fp);
         amass[ia] *= 1822.89;
         fread(&zv[ia],sizeof(REAL),1,fp);

         fread(&lmax0,sizeof(int),1,fp);
         if (lmax0!=lmax[ia])
            printf("Error reading psp %d, lmax0=%d <> lmax=%d\n",ia,lmax0,lmax[ia]);

         fread(rc[ia],sizeof(REAL),lmax0+1,fp);
         fread( &vnlnrm[ia*9],    sizeof(REAL),lmmax[ia],fp);
      }
      Parallel_ibcast(0,1,&icube0);
      Parallel_ibcast(0,1,&nfft0);
      Parallel_rbcast(0,1,&unit0);
      Parallel_ibcast(0,1,&lmax0);
      Parallel_rbcast(0,lmax0+1,rc[ia]);

      fread( &vl[ia*nfft3d],   sizeof(REAL),nfft3d,fp);
      fread( &vnl[ia*nfft3d*9],sizeof(REAL),lmmax[ia]*nfft3d,fp);
      if (Parallel_taskid==0)
         fclose(fp);
   }
   Parallel_rbcast(0,nkatm,amass);
   Parallel_rbcast(0,nkatm,zv);
   Parallel_rbcast(0,nkatm*9,vnlnrm);

   /* Setup Ewald summation */
   rs=pow( (3.00*omega/fourpi),(1.00/3.00));
   zz=0.0;
   z =0.0;
   for (ii=0; ii<nion; ++ii)
   {
     q = zv[katm[ii]];
     zz += q*q;
     z += q;
   }

   cewald = -0.50*zz*(alpha/rs+gssum(nfft,vg)/omega)
            -0.50*(z*z-zz)*rcut*rcut*pi/omega;

   /* lattice vectors in the cordination space */
   ii = 2*ncut+1;
   nsh = ii*ii*ii;
   rcell = (REAL *) malloc(3*nsh*sizeof(REAL));
   l=0;
   for (i3=-ncut; i3<=ncut; ++i3)
   for (i2=-ncut; i2<=ncut; ++i2)
   for (i1=-ncut; i1<=ncut; ++i1)
   {
      rcell[l]         = i1*unita[0]+i2*unita[3]+i3*unita[6];
      rcell[l+nsh]     = i1*unita[1]+i2*unita[4]+i3*unita[7];
      rcell[l+nsh+nsh] = i1*unita[2]+i2*unita[5]+i3*unita[8];
      ++l;
   }

   /* allocate cpsd work array */
   //cpsdlwork = 2*(nfft*nfft+10) + 40*nfft;
   cpsdlwork  = (6+ispin)*n2ft3d;
   cpsdlwork += nion*(5*nfft+3);
   cpsdlwork += 7*ne[0]*ne[0];
   cpsdwork = (REAL *) malloc(cpsdlwork*sizeof(REAL));


  /* geometrical center of the cluster and center of mass */
  for (i=0; i<3; ++i) gc[i] = 0.0;
  for (i=0; i<3; ++i) cm[i] = 0.0;
  am=0.0;
  for (ii=0; ii<nion; ++ii)
  {
     for (i=0; i<3; ++i) gc[i] += r2[i+3*ii];
     for (i=0; i<3; ++i) cm[i] += amass[katm[ii]]*r2[i+3*ii];
     am += amass[katm[ii]];
   }
  for (i=0; i<3; ++i) gc[i] /= ((REAL) nion);
  for (i=0; i<3; ++i) cm[i] /= am;

    
   /* summary of input data */
   if (oprint)
   {
      printf("          ================ input data ========================\n");

      printf("\n number of processors used: %d\n", Parallel_np());
      if (mapping==1) printf(" parallel mapping:             1d slab\n");
      if (mapping==2) printf(" parallel mapping:          2d hilbert\n");
      if (mapping==3) printf(" parallel mapping:           2d hcurve\n");
       printf("\n options:\n");
      if (imove)
         printf("      ionic motion   = yes\n");
      else
         printf("      ionic motion   = no\n");
      printf("      electron spin  = %d\n",ispin);

      printf("\n elements involved in the cluster:\n");
      for (ia=0; ia<nkatm; ++ia)
      {
         printf("   %2d : %s  mass no.: %6.1f   core charge %4.1f  lmax = %1d\n",
                ia+1,atom[ia],(amass[ia]/1822.89),zv[ia],lmax[ia]);
         printf("           cutoff = ");
         for (l=0; l<=lmax[ia]; ++l) printf(FMT8p3,rc[ia][l]);
         printf("\n");
      }


      printf("\n atomic composition:\n");
      for (ia=0; ia<nkatm; ++ia)
        printf("     %s : %d",atom[ia],natm[ia]);
      printf("\n");
      printf("\n initial position of ions:\n");
      for (ii=0; ii<nion; ++ii)
      {
          printf("     %4d %2s ( %11.5f %11.5f %11.5f )\n",
              ii+1,atom[katm[ii]],r2[3*ii],r2[1+3*ii],r2[2+3*ii]);
      }
      printf("        G.C. ( %11.5f %11.5f %11.5f )\n",gc[0],gc[1],gc[2]);
      printf("       C.O.M.( %11.5f %11.5f %11.5f )\n",cm[0],cm[1],cm[2]);

      printf("\n supercell:\n");
      if (icube==1) printf("     lattice: simple cubic");
      if (icube==2) printf("     lattice: face-centered cubic");
      if (icube==3) printf("     lattice: body-centered cubic");
      printf("   size="); printf(FMT8p3,unit);
      printf(" volume="); printf(FMT10p1,omega);
      printf("\n");
      printf("                 ");
      printf("a1 = < "); printf(FMT10,unita[0],unita[1],unita[2]); printf(" >\n");
      printf("                 ");
      printf("a2 = < "); printf(FMT10,unita[3],unita[4],unita[5]); printf(" >\n");
      printf("                 ");
      printf("a3 = < "); printf(FMT10,unita[6],unita[7],unita[8]); printf(" >\n");
      printf("     reciprocal: ");
      printf("b1 = < "); printf(FMT10,unitg[0],unitg[1],unitg[2]); printf(" >\n");
      printf("                 ");
      printf("b2 = < "); printf(FMT10,unitg[3],unitg[4],unitg[5]); printf(" >\n");
      printf("                 ");
      printf("b3 = < "); printf(FMT10,unitg[6],unitg[7],unitg[8]); printf(" >\n");
      printf("     cutoff=");printf(FMT8p3,ecut);
      printf("  fft=%4d x %4d x %4d (%8d waves)\n",nfft,nfft,nfft,nfft3d);
      printf("     Ewald summation: cut radius=");printf(FMT8p3,rcut);
      printf("  and %3d\n",ncut);
      printf("\n");

      printf("\n technical parameters:\n");
      printf("      iterations = %d (inner = %d, outer = %d)\n",inner*outer,inner,outer);
      printf("      time step  =");printf(FMT8p3,dt);
      printf("      fictitious mass=");printf(FMT10p1,fmass);
      printf("\n");
      printf("      tolerances =");
      printf(FMTE8p3,tole); printf(" (energy)");
      printf(FMTE12p3,tolc); printf(" (electron)");
      printf(FMTE12p3,tolr); printf(" (ion)\n");
      printf("\n\n");
      printf("          ================ iteration =========================\n");
      message(2);
   }

/*****************************|  begin iteration |****************************/
   cpu2 = current_second();
   done = 0;
   icount=0;
   for (i=0; i<10; ++i) E[i] = 0.0;
   deltae = 0.0;
   deltac = 0.0;
   deltar = 0.0;
   while (!done) 
   {
      ++icount;


      cpsd(inner,imove,nion,ispin,ne,nfft,nfft3d,n2ft3d,
           unitg,omega,dt,fmass,katm,amass,zv,nsh,rcell,cewald,rcut,
           masker,G,tg,vc,vg,lmmax,vnlnrm,vl,vnl,r1,r2,fion,c1,c2,cpsi,dn,
           hml,lmd,cpsdwork,
           E,&deltae,&deltac,&deltar);

      if (oprint)
      {
         printf("%8d",inner*icount);
         printf(FMTE20p10,E[0]);
         printf(FMTE13p5,deltae);
         printf(FMTE13p5,deltac);
         printf(FMTE13p5,deltar);
         printf("\n");
      }
      if ((fabs(deltae)< tole) && (deltac < tolc) && (deltar < tolr))
      {
         if (oprint) printf("          *** tolerances ok.  iteration terminated.\n");
         done = 1;
      }
      if (icount >= outer)
      {
         if (oprint) printf("          *** arrived at the maximum iteration.  iteration terminated.\n");
         done = 1;
      }
   }


/*****************************|  EPILOGUE  |****************************/
   if (oprint) message(3);
   cpu3 = current_second();

   /*  check orthnormality */
   n1[0] = 0; 
   n2[0] = ne[0]-1; 
   n1[1] = ne[0]; 
   n2[1] = ne[0]+ne[1]-1; 
   if (oprint) printf("\n orthonormality\n");
   for (ms=0; ms<ispin; ++ms)
   for (i=n1[ms]; i<=n2[ms]; ++i)
   for (j=i;      j<=n2[ms]; ++j)
   {
      //w = gcdotc(nfft,&c1[i*n2ft3d],&c1[j*n2ft3d]);
      w = d3db_cc_dot(&c1[i*n2ft3d],&c1[j*n2ft3d]);
      if (oprint) {printf("%3d %3d %3d",ms+1,i-n1[ms]+1,j-n1[ms]+1);printf(FMTE18p7,w);printf("\n");}
   }

   /* diagnalization of hamiltonian matrix */
   tch1 = 'V';
   tch2 = 'U';
   tch3 = 'N';
   for (ms=0; ms<ispin; ++ms)
   {
      esyev(&tch1,&tch2, &ne[ms], 
            &hml[ms*ne[0]*ne[0]], &ne[ms], 
           &eig[ms*ne[0]],
            work,&lwork,&info);
      eigsrt(&eig[ms*ne[0]],&hml[ms*ne[0]*ne[0]],ne[ms]);

      egemm(&tch3,&tch3,&n2ft3d,&ne[ms],&ne[ms],
               &rone,
               &c2[ms*ne[0]*n2ft3d], &n2ft3d,
               &hml[ms*ne[0]*ne[0]], &ne[ms],
               &rzero,
               &c1[ms*ne[0]*n2ft3d],&n2ft3d);
   }



   /* geometrical center of the cluster and center of mass */
   for (i=0; i<3; ++i) gc[i] = 0.0;
   for (i=0; i<3; ++i) cm[i] = 0.0;
   am=0.0;
   for (ii=0; ii<nion; ++ii)
   {
      for (i=0; i<3; ++i) gc[i] += r2[i+3*ii];
      for (i=0; i<3; ++i) cm[i] += amass[katm[ii]]*r2[i+3*ii];
      am += amass[katm[ii]];
   }
   for (i=0; i<3; ++i) gc[i] /= ((REAL) nion);
   for (i=0; i<3; ++i) cm[i] /= am;

   if (oprint)
   {
      printf("\n\n          ================ summary of results ================\n");
      printf("\n final position of ions:\n");
      for (ii=0; ii<nion; ++ii)
      {
          printf("     %4d %2s ( %11.5f %11.5f %11.5f )\n",
              ii+1,atom[katm[ii]],r2[3*ii],r2[1+3*ii],r2[2+3*ii]);
      }
      printf("        G.C. ( %11.5f %11.5f %11.5f )\n",gc[0],gc[1],gc[2]);
      printf("       C.O.M.( %11.5f %11.5f %11.5f )\n",cm[0],cm[1],cm[2]);

      printf("\n total     energy    :"); printf(FMTE19p10,E[0]); printf(" (");printf(FMTE15p5,E[0]/nion); printf("/ion)\n");
      printf(" total orbital energy:"); printf(FMTE19p10,E[1]); printf(" (");printf(FMTE15p5,E[1]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" hartree   energy    :"); printf(FMTE19p10,E[2]); printf(" (");printf(FMTE15p5,E[2]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" exc-corr  energy    :"); printf(FMTE19p10,E[3]); printf(" (");printf(FMTE15p5,E[3]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" ion-ion   energy    :"); printf(FMTE19p10,E[4]); printf(" (");printf(FMTE15p5,E[4]/nion); printf("/ion)\n");

      printf("\n KS kinetic energy   :"); printf(FMTE19p10,E[5]); printf(" (");printf(FMTE15p5,E[5]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" KS V_l energy       :"); printf(FMTE19p10,E[6]); printf(" (");printf(FMTE15p5,E[6]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" KS V_nl energy      :"); printf(FMTE19p10,E[7]); printf(" (");printf(FMTE15p5,E[7]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" KS V_Hart energy    :"); printf(FMTE19p10,E[8]); printf(" (");printf(FMTE15p5,E[8]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" KS V_xc energy      :"); printf(FMTE19p10,E[9]); printf(" (");printf(FMTE15p5,E[9]/(ne[0]+ne[1])); printf("/electron)\n");
      printf(" Virial Ratio <V>/<T>:"); printf(FMTE19p10,-99.0); printf("\n");

      printf("\n orbital energies:\n");
      nn = ne[0] - ne[1];
      ev = 27.2116;
      for (i=0; i<nn; ++i)
      {
         printf(FMTE18p7,eig[i]); printf(" ("); printf(FMT8p3,eig[i]*ev); printf("eV)\n");
      }
      for (i=0; i<ne[1]; ++i)
      {
         printf(FMTE18p7,eig[i+nn]); printf(" ("); printf(FMT8p3,eig[i+nn]*ev); printf("eV) ");
         printf(FMTE18p7,eig[i+(ispin-1)*ne[0]]); printf(" ("); printf(FMT8p3,eig[i+(ispin-1)*ne[0]]*ev); printf("eV)\n");
      }
   }

   /* writeout  ELCOUT */
   if (Parallel_taskid==0)
   {
      fp = fopen("ELCOUT","wb");
      fwrite(&icube,sizeof(int),1,fp);
      fwrite(&nfft,sizeof(int),1,fp);
      fwrite(&unit,sizeof(REAL),1,fp);
      fwrite(&ispin,sizeof(int),1,fp);
      fwrite(ne,sizeof(int),2,fp);
      nn = (ne[0]+ne[1])*n2ft3d;
   }
   //fwrite(c1,sizeof(REAL),nn,fp);
   for (i=0; i<(ne[0]+ne[1]); ++i)
      d3db_c_write(fp,&c1[i*n2ft3d],c2,cpsi);
   if (Parallel_taskid==0)
      fclose(fp);

   /*  writeout ionic structure */
   nn = 3*nion;
   ecopy(&nn,&rzero,&zero,r1,&one);
   if (Parallel_taskid==0)
   {
      fp = fopen("IONOUT","w");
      fprintf(fp,"%d\n",nkatm);
      for (ia=0; ia<nkatm; ++ia)
      {
         fprintf(fp,"\'%s \' %d\n",fname[ia],lmax[ia]);
      }
      fprintf(fp,"%d\n",nion);
      for (ii=0; ii<nion; ++ii)
      {
         fprintf(fp,"%d ",katm[ii]+1);
         fprintf(fp,FMT1,r2[0+ii*3]); fprintf(fp," ");
         fprintf(fp,FMT1,r2[1+ii*3]); fprintf(fp," ");
         fprintf(fp,FMT1,r2[2+ii*3]); fprintf(fp," ");
         fprintf(fp,FMT1,r1[0+ii*3]); fprintf(fp," ");
         fprintf(fp,FMT1,r1[1+ii*3]); fprintf(fp," ");
         fprintf(fp,FMT1,r1[2+ii*3]); fprintf(fp," ");
         fprintf(fp,"\n");
      }
      fclose(fp);
   }




   /* deallocate memory */
   free(rcell);
   free(r1);
   free(r2);
   free(fion);
   free(katm);
   for (ia=0; ia<nkatm; ++ia) 
   {
      free(fname[ia]);
      free(atom[ia]);
      free(rc[ia]);
   }
   free(fname);
   free(atom);
   free(rc);

   free(vnl);
   free(vl);
   free(vnlnrm);

   free(natm);
   free(lmax);
   free(masker);
   free(masker1);
   free(G);
   free(tg);
   free(vc);
   free(vg);

   free(cpsdwork);
   free(work);
   free(eig);
   free(lmd);
   free(hml);
   free(xce);
   free(xcp);
   free(dn);
   free(cpsi);
   free(c2);
   free(c1);

   d3db_end();
   Parallel_end();

   cpu4 = current_second();

   /* report consumed cputime */
   t1=cpu2-cpu1;
   t2=cpu3-cpu2;
   t3=cpu4-cpu3;
   t4=cpu4-cpu1;
   av=t2/((REAL) (icount*inner));
   if (oprint)
   {
      printf("\n -----------------\n");
      printf(" cputime in seconds\n");
      printf(" prologue    : "); printf(FMTE15p5,t1);printf("\n");
      printf(" main loop   : "); printf(FMTE15p5,t2);printf("\n");
      printf(" epilogue    : "); printf(FMTE15p5,t3);printf("\n");
      printf(" total       : "); printf(FMTE15p5,t4);printf("\n");
      printf(" cputime/step: "); printf(FMTE15p5,av);printf("\n");
      message(4);
   }

}
