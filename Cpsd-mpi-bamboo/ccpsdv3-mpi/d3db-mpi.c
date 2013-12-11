/*
 * d3db-mpi.c
 *
 *  Created on: Dec 5, 2013
 *      Author: bylaska
 *
 *
, i i i i i t i i ; , ; ; ; , ; ; i i ; ; : i ; , ; : , , : , , , , , : : : , : : ; ; ; i ; i , , :
; i i ; i ; i i , ; , ; , i , t i , f E ; , i ; ; . , , , , : : : : , , : : : : : : : : , , : , , .
, ; ; : , ; : ; , , , i , , i i i ; j i j L f i , i ; i : t , : : : . , : : : , : : : : , : : : : .
, ; j i i ; i ; i ; t , , G D t i i t i t t f D : # G L f L D G , ; , , . : : : : : : : : : : : : .
i t i t t t ; i i i i G ; G , ; E t i i i L t i K G # W E L D D L K : , , : : : : : . . : . : : : .
i i i i t i i ; t i ; i t j t D # W K t j j i j L E # G W G # D i L : : : : : : : : : : . : : : : :
; j i t i ; i ; ; ; t j i K # L f t i j t i j i ; E # # # W G L # ; i , , ; , , ; : ; : . : . , : :
; ; ; , i ; i t t i j , t t W j j j G i ; D , ; j ; j t # # E W W ; i i , , ; , ; , : , : : , , , ,
i ; i t t j L f L t t f # K j K i t t t ; i ; : ; j , D G D W K # L W , . , , : , : : : , , , : : :
; t j j j t f j f f j i f j i t , i ; ; i i G ; f i t f f W L E f j t : . : : : : : : , : , : : : :
f j L f j t i t i t ; , : j ; t K t L j i t # # ; ; , K L , . D f , : : . , : , : : : . . : . : : ,
i i ; t i i ; ; i i ; , ; ; ; W K f K G i W t ; ; ; t G i t j j , t D i , , ; : : . : . . : : . . .
; ; i ; ; i ; i i i t ; ; ; i i K W W K # i ; , i , , L G E G K L E i j G D t i : : . . . . : : : .
; i , ; i i i ; i t ; ; ; ; ; i t # j t ; j ; ; ; i f D ; W D D K W E K D W i : : : : : : : , : , :
i i ; ; i t i t i i i ; ; ; i ; t E t i t ; ; i i ; W D # # # K E E E j ; : t : , : : : : : : . : ;
; t t ; i t ; t ; i ; ; ; ; ; ; W j ; ; ; ; i ; i : f D # # W D K f f , ; , i j , : : : : : , : : :
i ; t i i i t t ; i i ; i ; ; t # i i ; ; i ; , ; E j E # # E D D t j L L D j f : : , , , : , : ; :
t i i i i t i i i i i i i i i i # f j ; ; ; t i : j t ; E E j i i f i f G t j ; , , : : : : : : : :
i i i t t t i i t t i i i i ; i # t , t f K K D L j L i E G D G f f t , ; ; ; ; , : : : : , : : : ,
i i t t t t t t t t j j t t i j # W K # i t t , f i t t t j t i i ; , ; ; ; ; ; ; , ; , , , , , , ,
t t i i t t i j t t t t t t i # W i ; ; i t i t i i i i i i ; i ; i i i i i ; ; , ; , ; , i ; ; , ;
t t t t i t t j j t i t t t i # i i i t i i i i i i i i i t i ; i ; ; ; ; ; , ; i ; ; ; ; ; i ; , ,
j j t j t t j t t t t j t j W L i ; i i i i t t i i i t i i i ; i i i ; ; ; ; i t i ; t , i i ; ; ;
j j j j t t t t t t t i W W # t i ; i i i i i i i i i ; i i i i i ; , ; i i i ; i f L j j D t ; , ;
j j t t t t t t t ; i f K # ; i j i i i j t j L W W K K D ; f ; i ; i , t D G L W L j G # # W D K D
t j j L f t t t t t i t # i G # K K # # G W W K W E W W E K E K K W D E E E E K K K L # E D E K W W
W D f K K D i t t f i W # # W W # W # # # K W W K K E K E K E K E K E E K K E D t W D D L L G G D K
# W # # # # # # # # D # W K W W K K W W W W W K W K E E K K K E K E K K K K # E f K # W W D E E W #
# W # # # W W W # K # # W K K # W W W # # # # # # # # # W K # # # E W # # # # K W W W # # # # # # #
# # # # # # # # # # # # # # # K # # # # # # # # # # # # # # # # W K W K E K W W E K K K K W G K K E
W W W W # # W K W # # # W K E E K E K K K K K E G D K D D D D E D E E D D E K E E D D D G D E D K K
W K W W K K E W K W L W # # # # # # # K K W W # W D E G G L E E L D D E D E E E E E D D G D D G K K
j t L f E G f t W # W K # W W # W E ; W # K K W G E E K D W W W D K K D K G G D K E D E D K D D D K
G E G G f D L D L L f L L L D E E E K K K W K W K W G K K E f D G G E D D D G D E D f D f G K D G D
G E E D D D E E D D D D G D D D L G D D K E E E K K ; D G D D D i D f E L D D f L G G G L L D G D D
*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "float.h"
#include "Parallel.h"
#include "d3db.h"
#include "hilbert.h"
#include "hcurve.h"

/********************/
/* static variables */
/********************/
static int nq,nx,ny,nz;
static int nfft3d,n2ft3d,nfft3d_map,n2ft3d_map;
static int zplane_size;
static int mapping, mapping2d;

/* slab mapping */
static int *q_map, *p_map, *k_map;
static int *i1_start, *i2_start, *iq_to_i1, *iq_to_i2;

/* hilbert/h-curve mapping */
static int *q_map1,*p_map1,nq1;
static int *q_map2,*p_map2,nq2;
static int *q_map3,*p_map3,nq3;
static int *h_i1_start[6],*h_i2_start[6],*h_iq_to_i1[6],*h_iq_to_i2[6];

/* fft workspace variables */
static REAL *tmpx,*tmpy,*tmpz;


/* data fetch routines */
int d3db_nfft3d() {return nfft3d;}
int d3db_nfft3d_map() {return nfft3d_map;}
int d3db_n2ft3d() {return n2ft3d;}
int d3db_n2ft3d_map() {return n2ft3d_map;}
int d3db_nq() {return nq; }
int d3db_nx() {return nx;}
int d3db_ny() {return ny;}
int d3db_nz() {return nz;}
int d3db_zplane_size() {return zplane_size;}


/******************************************
 *                                        *
 *          generate_map_indexes          *
 *                                        *
 ******************************************
*/
static void generate_map_indexes(int taskid, int np, int ny, int nz,
		                 int p_map[], int q_map[], int *nq_out)
{
	int i,j,k,p,q,nq,nq1,nq2,rmdr1;
	int *indx_proc,*indx_q,*tmp_p;

	indx_proc = (int *) malloc(ny*nz*sizeof(int));
	indx_q    = (int *) malloc(ny*nz*sizeof(int));
	tmp_p     = (int *) malloc(ny*nz*sizeof(int));

	for (i=0; i<(ny*nz); ++i)
		tmp_p[i] = p_map[i];
	nq1 = (ny*nz)/np;
	rmdr1 = (ny*nz)%np;
	nq2 = nq1;
	if (rmdr1>0) ++nq2;
	nq = 0;
	p = 0;
	q = 0;
	for (i=0; i<(ny*nz); ++i)
	{
		indx_proc[i] = p;
		indx_q[i] = q;
		if (taskid==p) ++nq;
		++q;
		if (q>=nq2)
		{
		   q = 0;
		   ++p;
		   p = p%np;
		   if (p>=rmdr1) nq2 = nq1;
		}
	}

	for (k=0; k<nz; ++k)
	for (j=0; j<ny; ++j)
	{
	   i = tmp_p[j+k*ny];
	   p = indx_proc[i];
	   q = indx_q[i];
	   p_map[j+k*ny] = p;
	   q_map[j+k*ny] = q;
	}

	free(tmp_p);
	free(indx_q);
	free(indx_proc);

	*nq_out = nq;
}

/*******************************************
 *                                         *
 *              mapping_init               *
 *                                         *
 *******************************************
*/
static void mapping_init()
{
	int np     = Parallel_np();
	int taskid = Parallel_taskid();

	/* slab mapping */
	if (mapping==1)
	{
	   int p,q,k;
	   q_map = (int *) malloc(nz*sizeof(int));
	   p_map = (int *) malloc(nz*sizeof(int));
	   k_map = (int *) malloc(nz*sizeof(int));

	   /* cyclic */
	   p = 0; q = 0;
	   for (k=0; k<nz; ++k)
	   {
	      q_map[k] = q;
	      p_map[k] = p;
	      if (p==taskid) nq = q+1;
          ++p;
          if (p>=np){p = 0; ++q;}
	   }
	   for (k=0; k<nz; ++k)
	      if (p_map[k]==taskid)
		     k_map[q_map[k]] = k;
	   nfft3d = (nx/2+1)*ny*nq;
           n2ft3d = 2*nfft3d;
	   nfft3d_map = nfft3d;
	   n2ft3d_map = n2ft3d;
	}

	/* Hilbert mapping */
	else
	{
	   q_map1 = (int *) malloc(ny*nz*sizeof(int));
	   p_map1 = (int *) malloc(ny*nz*sizeof(int));
	   q_map2 = (int *) malloc(nz*(nx/2+1)*sizeof(int));
	   p_map2 = (int *) malloc(nz*(nx/2+1)*sizeof(int));
	   q_map3 = (int *) malloc(ny*(nx/2+1)*sizeof(int));
	   p_map3 = (int *) malloc(ny*(nx/2+1)*sizeof(int));

           if (mapping2d==1)
           {
	      hilbert2d_map_c(ny,nz,    p_map1);
	      hilbert2d_map_c(nz,nx/2+1,p_map2);
	      hilbert2d_map_c(nx/2+1,ny,p_map3);
           }
           else
           {
              hcurve_map_c(ny,nz,    p_map1);
	      hcurve_map_c(nz,nx/2+1,p_map2);
	      hcurve_map_c(nx/2+1,ny,p_map3);
           }
	   generate_map_indexes(taskid,np,ny,nz,    p_map1,q_map1,&nq1);
	   generate_map_indexes(taskid,np,nz,nx/2+1,p_map2,q_map2,&nq2);
	   generate_map_indexes(taskid,np,nx/2+1,ny,p_map3,q_map3,&nq3);

		nfft3d = (nx/2+1)*nq1;
		if ((ny*nq2)>nfft3d) nfft3d = ny*nq2;
		if ((nz*nq3)>nfft3d) nfft3d = nz*nq3;
		n2ft3d = 2*nfft3d;

		nfft3d_map = nz*nq3;
		n2ft3d_map = (nx+2)*nq1;
	}
}


static void d3db_fftbx_sub(int n, int nx, int nxh,
		                   REAL tmpx[], REAL A[])
{
   int i;
   REAL *atmp;

   atmp = A;
   for (i=0; i<n; ++i)
   {
	   erfftb(&nx,atmp,tmpx);
	   atmp = &atmp[2*nxh];
   }
}

static void d3db_fftby_sub2(int n, int ny,
                            REAL tmpy[], REAL A[])
{
	int i;
	REAL *atmp;
	atmp = A;
	for (i=0; i<n; ++i)
	{
		ecfftb(&ny,atmp,tmpy);
		atmp = &atmp[2*ny];
	}
}

static void d3db_fftbz_sub2(int n, int nz,
                            REAL tmpz[], REAL A[])
{
	int i;
	REAL *atmp;

	atmp = A;
	for (i=0; i<n; ++i)
	{
		ecfftb(&nz,atmp,tmpz);
		atmp = &atmp[2*nz];
	}
}


static void cshift1_fftb(int nx, int ny, int nq, int ne, REAL A[])
{
   int i,j,indx;

   indx = 0;
   for (j=0; j<(ny*nq*ne); ++j)
   {
	   for (i=1; i<nx; ++i)
          A[indx+i] = A[indx+i+1];
      indx = indx + (nx+2);
   }
}


static void zeroend_fftb(int nx, int ny,int nq, int ne, REAL A[])
{
   int i,indx;

   indx  = nx;
   for (i=0; i<(ny*nq*ne); ++i)
   {
      A[indx]   = 0.0;
      A[indx+1] = 0.0;
      indx      = indx + (nx+2);
   }
}

static void cshift_fftf(int nx, int ny, int nq, int ne, REAL A[])
{
   int i,j,indx;

   indx = 0;
   for (j=0; j<(ny*nq*ne); ++j)
   {
      for (i=nx; i>=2; --i)
	     A[indx+i] = A[indx+i-1];
	  A[indx+1]    = 0.0;
      A[indx+nx+1] = 0.0;
	  indx += (nx+2);
   }
}



/***********************************
 *					               *
 *	          d3db_qtok   		   *
 *					               *
 ***********************************
*/
void d3db_qtok(int q, int *k) {*k = k_map[q];}


/***********************************
 *					               *
 *	          d3db_ktoqp  		   *
 *					               *
 ***********************************
*/
void d3db_ktoqp(int k, int *q, int *p) {*q = q_map[k]; *p = p_map[k];}

/***********************************
 *					               *
 *	       d3db_ijktoindexp	       *
 *					               *
 ***********************************
*/
void d3db_ijktoindexp(int i, int j, int k,
		              int *indx, int *p)
{
   int q;

   /**** slab mapping ***/
   if (mapping==1)
   {
      q  = q_map[k];
      *p = p_map[k];
      *indx = i + j*(nx/2+1) + q*(nx/2+1)*ny;
   }

   /**** Hilbert mapping ****/
   else
   {
      q  = q_map3[i+j*(nx/2+1)];
      *p = p_map3[i+j*(nx/2+1)];
      *indx = k + q*nz;
   }
}


/***********************************
 *                                 *
 *           d3db_ijktoindex1p     *
 *                                 *
 ***********************************
*/
void d3db_ijktoindex1p(int i, int j, int k,
		               int *indx, int *p)
{
   int q;

   /**** slab mapping ***/
   if (mapping==1)
   {
      q  = q_map[j];
      *p = p_map[j];
      *indx = i + k*(nx/2+1) + q*(nx/2+1)*nz;
   }
   /**** hilbert mapping ****/
   else
   {
      q = q_map2[k+i*nz];
      *p = p_map2[k+i*nz];
      *indx = j + q*ny;
   }
}

/***********************************
 *                                 *
 *       d3db_ijktoindex2p         *
 *                                 *
 ***********************************
*/
void d3db_ijktoindex2p(int i,int j,int k,int *indx,int *p)
{
   int q;

   /* slab mapping */
   if (mapping==1)
   {
      q = q_map[j];
      *p = p_map[j];
      *indx = i + k*(nx+2) + q*(nx+2)*ny;
   }

   /* Hilbert mapping */
   else
   {
      q = q_map1[j+k*ny];
      *p = p_map1[j+k*ny];
      *indx = i + q*(nx+2);
   }
}



/*****************************************
 *                                       *
 *        d3db_c_transpose_jk_init       *
 *                                       *
 *****************************************
*/
static void d3db_c_transpose_jk_init()
{
	int index1,index2,it,proc_to,proc_from,qhere,phere,qto,pto,qfrom,pfrom,itmp,i,j,k;

	int taskid = Parallel_taskid();
	int np     = Parallel_np();

	iq_to_i1 = (int *) malloc((nx/2+1)*ny*nq*sizeof(int));
	iq_to_i2 = (int *) malloc((nx/2+1)*ny*nq*sizeof(int));
	i1_start = (int *) malloc((nz+1)*sizeof(int));
	i2_start = (int *) malloc((nz+1)*sizeof(int));

    index1 = 0;
    index2 = 0;
    for (it=0; it<np; ++it)
    {
       proc_to   = (taskid+it)%np;
       proc_from = (taskid-it+np)%np;
       i1_start[it] = index1;
       i2_start[it] = index2;

       for (k=0; k<nz; ++k)
       for (j=0; j<ny; ++j)
       {
          /**** packing scheme ****/
          d3db_ktoqp(k,&qhere,&phere);
          d3db_ktoqp(j,&qto,&pto);
          if ((phere==taskid) && (pto==proc_to))
          {
        	 for (i=0; i<(nx/2+1); ++i)
        	 {
                itmp = i + j*(nx/2+1) + qhere*(nx/2+1)*ny;
                iq_to_i1[itmp] = index1;
                ++index1;
        	 }
          }

          /**** unpacking scheme ****/
          d3db_ktoqp(j,&qhere,&phere);
          d3db_ktoqp(k,&qfrom,&pfrom);
          if ((phere==taskid) && (pfrom==proc_from))
          {
        	 for (i=0; i<(nx/2+1); ++i)
        	 {
                itmp = i + k*(nx/2+1) + qhere*(nx/2+1)*ny;
                iq_to_i2[itmp] = index2;
                ++index2;
        	 }
          }
       }
    }
    i1_start[np] = index1;
    i2_start[np] = index2;
}


/*****************************************
 *                                       *
 *       d3db_c_transpose_ijk_init       *
 *                                       *
 *****************************************
*/
static void d3db_c_transpose_ijk_init()
{

	int itmp,index1,index2,it,proc_to,proc_from,pto,qto,phere,qhere,i,j,k;

	int taskid = Parallel_taskid();
	int np     = Parallel_np();

   /*********************************************************
    **** map1to2 mapping - done - transpose operation #1 ****
    *********************************************************/
   h_iq_to_i1[0] = (int *) malloc((nx/2+1)*nq1*sizeof(int));
   h_iq_to_i2[0] = (int *) malloc(ny*nq2*sizeof(int));
   h_i1_start[0] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[0] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[0][it] = index1;
      h_i2_start[0][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {
         /**** packing scheme ****/
         phere = p_map1[j+k*ny];
         qhere = q_map1[j+k*ny];
         pto   = p_map2[k+i*nz];
         qto   = q_map2[k+i*nz];

         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = i + qhere*(nx/2+1);
            h_iq_to_i1[0][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = j + qto*ny;
            h_iq_to_i2[0][itmp-1] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[0][np] = index1;
   h_i2_start[0][np] = index2;


   /*********************************************************
    **** map2to3 mapping - done - transpose operation #2 ****
    *********************************************************/
   h_iq_to_i1[1] = (int *) malloc(ny*nq2*sizeof(int));
   h_iq_to_i2[1] = (int *) malloc(nz*nq3*sizeof(int));
   h_i1_start[1] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[1] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[1][it] = index1;
      h_i2_start[1][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {
         /**** packing scheme ****/
         phere = p_map2[k+i*nz];
         qhere = q_map2[k+i*nz];
         pto   = p_map3[i+j*(nx/2+1)];
         qto   = q_map3[i+j*(nx/2+1)];

         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = j + qhere*ny;
            h_iq_to_i1[1][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = k + qto*nz;
            h_iq_to_i2[1][itmp] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[1][np] = index1;
   h_i2_start[1][np] = index2;


   /*********************************************************
    **** map3to2 mapping - done - transpose operation #3 ****
    *********************************************************/
   h_iq_to_i1[2] = (int *) malloc(nz*nq3*sizeof(int));
   h_iq_to_i2[2] = (int *) malloc(ny*nq2*sizeof(int));
   h_i1_start[2] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[2] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[2][it] = index1;
      h_i2_start[2][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {
         /**** packing scheme ****/
         phere = p_map3[i+j*(nx/2+1)];
         qhere = q_map3[i+j*(nx/2+1)];
         pto   = p_map2[k+i*nz];
         qto   = q_map2[k+i*nz];

         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = k + qhere*nz;
            h_iq_to_i1[2][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = j + qto*ny;
            h_iq_to_i2[2][itmp] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[2][np] = index1;
   h_i2_start[2][np] = index2;


   /*********************************************************
    **** map2to1 mapping - done - transpose operation #4 ****
    *********************************************************/
   h_iq_to_i1[3] = (int *) malloc(ny*nq2*sizeof(int));
   h_iq_to_i2[3] = (int *) malloc((nx/2+1)*nq1*sizeof(int));
   h_i1_start[3] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[3] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[3][it] = index1;
      h_i2_start[3][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {
         /**** packing scheme ****/
         phere = p_map2[k+i*nz];
         qhere = q_map2[k+i*nz];
         pto   = p_map1[j+k*ny];
         qto   = q_map1[j+k*ny];

         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = j + qhere*ny;
            h_iq_to_i1[3][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = i + qto*(nx/2+1);
            h_iq_to_i2[3][itmp] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[3][np] = index1;
   h_i2_start[3][np] = index2;


   /***********************************************************
    **** map1to3 mapping  - done - transpose operation # 5 ****
    ***********************************************************/
   h_iq_to_i1[4] = (int *) malloc((nx/2+1)*nq1*sizeof(int));
   h_iq_to_i2[4] = (int *) malloc(nz*nq3*sizeof(int));
   h_i1_start[4] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[4] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[4][it] = index1;
      h_i2_start[4][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {
         /**** packing scheme ****/
         phere = p_map1[j+k*ny];
         qhere = q_map1[j+k*ny];
         pto   = p_map3[i+j*(nx/2+1)];
         qto   = q_map3[i+j*(nx/2+1)];

         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = i + qhere*(nx/2+1);
            h_iq_to_i1[4][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = k + qto*nz;
            h_iq_to_i2[4][itmp] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[4][np] = index1;
   h_i2_start[4][np] = index2;


   /*************************
    **** map3to1 mapping ****
    *************************/
   h_iq_to_i1[5] = (int *) malloc(nz*nq3*sizeof(int));
   h_iq_to_i2[5] = (int *) malloc((nx/2+1)*nq1*sizeof(int));
   h_i1_start[5] = (int *) malloc((np+1)*sizeof(int));
   h_i2_start[5] = (int *) malloc((np+1)*sizeof(int));

   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      h_i1_start[5][it] = index1;
      h_i2_start[5][it] = index2;

      for (k=0; k<nz; ++k)
      for (j=0; j<ny; ++j)
      for (i=0; i<(nx/2+1); ++i)
      {

         /**** packing scheme ****/
         phere = p_map3[i+j*(nx/2+1)];
         qhere = q_map3[i+j*(nx/2+1)];
         pto   = p_map1[j+k*ny];
         qto   = q_map1[j+k*ny];


         if ((phere==taskid) && (pto==proc_to))
         {
            itmp = k + qhere*nz;
            h_iq_to_i1[5][itmp] = index1;
            ++index1;
         }

         /**** unpacking scheme ****/
         if ((pto==taskid) && (phere==proc_from))
         {
            itmp = i + qto*(nx/2+1);
            h_iq_to_i2[5][itmp] = index2;
            ++index2;
         }
      }
   }
   h_i1_start[5][np] = index1;
   h_i2_start[5][np] = index2;

}


/*****************************************
 *                                       *
 *             d3db_fft_init             *
 *                                       *
 *****************************************/
static void d3db_fft_init()
{
   tmpx = (REAL *) malloc(2*(2*nx+15)*sizeof(REAL));
   tmpy = (REAL *) malloc(2*(2*ny+15)*sizeof(REAL));
   tmpz = (REAL *) malloc(2*(2*nz+15)*sizeof(REAL));
   erffti(&nx,tmpx);
   ecffti(&ny,tmpy);
   ecffti(&nz,tmpz);
}

/*****************************************
 *                                       *
 *              d3db_fft_end             *
 *                                       *
 *****************************************/
static void d3db_fft_end()
{
	free(tmpx);
	free(tmpy);
	free(tmpz);
}



/***********************************
 *                                 *
 *       d3db_c_transpose_jk       *
 *                                 *
 ***********************************

      This routine performs the operation
               A(i,k,j) <- A(i,j,k)

      np = the number of worker nodes
      proc#=0...(np-1)

       this transpose uses more buffer space
       then transpose2
*/

static void d3db_c_transpose_jk(REAL A[], REAL tmp1[], REAL tmp2[])
{
   int i,it,proc_to,proc_from,msgtype,msglen;
   int reqcnt;
   MPI_Request *request;

   int taskid = Parallel_taskid();
   int np     =  Parallel_np();

   request = (MPI_Request *) malloc(np*sizeof(MPI_Request));

   /* pack A(i) array */
   for (i=0; i<nfft3d; ++i)
   {
      tmp1[2*iq_to_i1[i]]   = A[2*i];
      tmp1[2*iq_to_i1[i]+1] = A[2*i+1];
   }

   /* it = 0, transpose data on same thread */
   msglen = i2_start[1] - i2_start[0];
   for (i=0; i<(2*msglen); ++i)
	   tmp2[2*i2_start[0]+i] = tmp1[2*i1_start[0]+i];


   /* receive packed array data */
   reqcnt = 0;
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;

      msgtype = 9;
      msglen = i2_start[it+1]-i2_start[it];

      if (msglen>0)
      {
         if (MPI_Irecv(&tmp2[2*i2_start[it]],2*msglen,MPI_REAL_PRECISION,
                        proc_from,msgtype,
                        MPI_COMM_WORLD,&request[reqcnt]) != MPI_SUCCESS)
            printf("d3db_c_transpose_jk error: MPI_Irecv failed\n");
         ++reqcnt;
      }
   }

   /* send packed array to other processors */
   for (it=1; it<np; ++it)
   {
      /* synchronous send of tmp */
      proc_to = (taskid+it)%np;
      msglen  = i1_start[it+1] - i1_start[it];
      msgtype   = 9;

      if (msglen>0)
         if (MPI_Send(&tmp1[2*i1_start[it]], 2*msglen,MPI_REAL_PRECISION,
                      proc_to,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
            printf("d3db_c_transpose_jk error: MPI_Send failed\n");
   }

   /* wait for completion of mp_send, also do a sync */
   if (np>1)
      if (MPI_Waitall(reqcnt,request,MPI_STATUS_IGNORE) != MPI_SUCCESS)
         printf("d3db_c_transpose_jk error: MPI_Waitall failed\n");

   /* unpack A(i) array */
   for (i=0; i<nfft3d; ++i)
   {
      A[2*i]   = tmp2[2*iq_to_i2[i]];
      A[2*i+1] = tmp2[2*iq_to_i2[i]+1];
   }
   free(request);
}



/***********************************
 *					               *
 *	   d3db_c_transpose_ijk		   *
 *					               *
 ***********************************

This routine performs the operation
        A(i,k,j) <- A(i,j,k)

np = the number of worker nodes
proc#=0...(np-1)

This transpose uses more buffer space then transpose2

*/
static void d3db_c_transpose_ijk(int op, REAL A[], REAL tmp1[], REAL tmp2[])
{

   int i,it,proc_to,proc_from,msgtype,msglen,nnfft3d;
   int reqcnt;
   MPI_Request *request;
   MPI_Status  *status;

   int taskid = Parallel_taskid();
   int np     =  Parallel_np();

   request = (MPI_Request *) malloc(np*sizeof(MPI_Request));
   status  = (MPI_Status *)  malloc(4*np*sizeof(MPI_Status));


   /* pack A(i) array */
   if ((op==1)||(op==5)) nnfft3d = (nx/2+1)*nq1;
   if ((op==2)||(op==4)) nnfft3d = (ny)    *nq2;
   if ((op==3)||(op==6)) nnfft3d = (nz)    *nq3;

   for (i=0; i<nnfft3d; ++i)
   {
      tmp1[2*h_iq_to_i1[op][i]]   = A[2*i];
      tmp1[2*h_iq_to_i1[op][i]+1] = A[2*i+1];
   }

   /* it = 0, transpose data on same thread */
   msglen = h_i2_start[op][1] - h_i2_start[op][0];
   for (i=0; i<(2*msglen); ++i)
      tmp2[2*h_i2_start[op][0]+i] = tmp1[2*h_i1_start[op][0]+i];

   /* receive packed array data */
   reqcnt = 0;
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msgtype = 29;
      msglen = h_i2_start[op][it+1] = h_i2_start[op][it];

      if (msglen>0)
      {
         if (MPI_Irecv(&tmp2[2*h_i2_start[op][it]],2*msglen,MPI_REAL_PRECISION,
        		       proc_from,msgtype,
        		       MPI_COMM_WORLD,&request[reqcnt])!=MPI_SUCCESS)
            printf("d3db_c_transpose_ijk error: MPI_Irecv failed\n");
         ++reqcnt;
      }
   }

   /* send packed array to other processors */
   for (it=1; it<np; ++it)
   {
      /* synchronous send of tmp */
         proc_to = (taskid+it)%np;
         msglen  = h_i1_start[op][it+1] - h_i1_start[op][it];
         msgtype = 29;

         if (msglen>0)
            if (MPI_Send(&tmp1[2*h_i1_start[op][it]],2*msglen,MPI_REAL_PRECISION,
                         proc_to,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
               printf("d3db_c_transpose_ijk error: MPI_Send failed\n");
   }

   /* wait for completion of mp_send, also do a sync */
   if (np>1)
      if (MPI_Waitall(reqcnt,request,status) != MPI_SUCCESS)
         printf("d3db_c_transpose_ijk error: MPI_Waitall failed\n");

   /* unpack A(i) array */
   if ((op==4)||(op==6)) nnfft3d = (nx/2+1)*nq1;
   if ((op==1)||(op==3)) nnfft3d = (ny)    *nq2;
   if ((op==2)||(op==5)) nnfft3d = (nz)    *nq3;

   for (i=0; i<nnfft3d; ++i)
   {
      A[2*i]   = tmp2[2*h_iq_to_i2[op][i]];
      A[2*i+1] = tmp2[2*h_iq_to_i2[op][i]+1];
   }
   free(request);
   free(status);
}


/***********************************
 *                                 *
 *      d3db_t_transpose_ijk       *
 *                                 *
 ***********************************
This routine performs the operation
         A(i,k,j) <- A(i,j,k)

np = the number of worker nodes
proc#=0...(np-1)

this transpose uses more buffer space then transpose2
*/
static void d3db_t_transpose_ijk(int op, REAL A[], REAL tmp1[], REAL tmp2[])
{
   int i,it,proc_to,proc_from,msgtype,msglen,nnfft3d;
   int reqcnt;
   MPI_Request *request;
   MPI_Status  *status;

   int taskid = Parallel_taskid();
   int np     =  Parallel_np();

   request = (MPI_Request *) malloc(np*sizeof(MPI_Request));
   status  = (MPI_Status *)  malloc(4*np*sizeof(MPI_Status));

   /* pack A(i) array */
   if ((op==1)||(op==5)) nnfft3d = (nx/2+1)*nq1;
   if ((op==2)||(op==4)) nnfft3d = (ny)    *nq2;
   if ((op==3)||(op==6)) nnfft3d = (nz)    *nq3;

   for (i=0; i<nnfft3d; ++i)
      tmp1[h_iq_to_i1[op][i]] = A[i];

   /* it = 0, transpose data on same thread */
   msglen = h_i2_start[op][1] - h_i2_start[op][0];
   for (i=0; i<msglen; ++i)
      tmp2[h_i2_start[op][0]+i] = tmp1[h_i1_start[op][0]+i];

   /* receive packed array data */
   reqcnt = 0;
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msgtype = 39;
      msglen = h_i2_start[op][it+1] + h_i2_start[op][it];

      if (msglen>0)
      {
         if (MPI_Irecv(&tmp2[h_i2_start[op][it]],msglen,MPI_REAL_PRECISION,
        		       proc_from,msgtype,
        		       MPI_COMM_WORLD,&request[reqcnt])!=MPI_SUCCESS)
        	 printf("d3db_t_transpose_ijk error: MPI_Irecv failed\n");
         ++reqcnt;
      }
   }

   /* send packed array to other processors */
   for (it=1; it<np; ++it)
   {
      /* synchronous send of tmp */
      proc_to = (taskid+it)%np;
      msglen  = h_i1_start[op][1] - h_i1_start[op][0];
      msgtype = 39;

      if (msglen>0)
         if (MPI_Send(&tmp1[h_i1_start[op][it]],msglen,MPI_REAL_PRECISION,
                      proc_to,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
            printf("d3db_t_transpose_ijk error: MPI_Send failed\n");
   }

   /* wait for completion of mp_send, also do a sync */
   if (np>1)
      if (MPI_Waitall(reqcnt,request,status) != MPI_SUCCESS)
         printf("d3db_t_transpose_ijk error: MPI_Waitall failed\n");

   /* unpack A(i) array */
   if ((op==4)||(op==6)) nnfft3d = (nx/2+1)*nq1;
   if ((op==1)||(op==3)) nnfft3d = (ny)    *nq2;
   if ((op==2)||(op==5)) nnfft3d = (nz)    *nq3;

   for (i=0; i<nnfft3d; ++i)
	   A[i] = tmp2[h_iq_to_i2[op][i]];
   free(request);
   free(status);
}





/* |------||------||------||------|\\//|------||------||------||------||------|
    //\\//\\//\\//\\//\\///\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\/
    \\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
    //\\//\\//\\//\\//\\///\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\/
    \\//\\//\\//\\//\\/(    )\\//\\//\\//\\//\\//\\//\\//\\//\\//\\///\\//\\//
    //\\//\\//\\//\\//\\(oo)//\\//\\//\\//\\//\\//\\//\\//\\//\//\\//\\//\\//\
    \\//\\//\\/)\.-----/(O O)\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\
    //\\//\\//# ;       / u//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\\//\\//
    \\//\\//\\//(  .   |} )\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\/
    //\\//\\//\\/|/ `.;|/;//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
    \\//\\//\\//\"     " "\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    \\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
    //\\//\\//\\//\\//\\///\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\/
   |------||------||------||------|\\//|------||------||------||------||------|*/





/*****************************************
 *                                       *
 *              d3db_end                 *
 *                                       *
 *****************************************/
void d3db_end()
{
   int i;
   //d3db_timereverse_end();
   d3db_fft_end();

   if (mapping==1)
   {
      free(q_map);
      free(p_map);
      free(k_map);

      free(i1_start);
      free(i2_start);
      free(iq_to_i1);
      free(iq_to_i2);
   }
   else
   {
      free(q_map1);
      free(p_map1);
      free(q_map2);
      free(p_map2);
      free(q_map3);
      free(p_map3);
      for (i=0; i<6; ++i)
      {
         free(h_i1_start[i]);
         free(h_i2_start[i]);
         free(h_iq_to_i1[i]);
         free(h_iq_to_i2[i]);
      }
   }
}

/*****************************************
 *                                       *
 *               d3db_init               *
 *                                       *
 *****************************************

*/
void d3db_init(int nx_in, int ny_in, int nz_in, int map_in)
{
   int np     = Parallel_np();
   int taskid = Parallel_taskid();

   nx = nx_in;
   ny = ny_in;
   nz = nz_in;
   mapping = map_in;
   mapping2d = 1;
   if (mapping==3)
   {
      mapping = 2;
      mapping2d=2;
   }
   mapping_init();
   if (mapping==1) d3db_c_transpose_jk_init();
   if (mapping==2) d3db_c_transpose_ijk_init();

   //d3db_c_timereverse_init();
   d3db_fft_init();
}


/***********************************
 *			                       *
 *	     d3db_cr_fft3b       	   *
 *			                       *
 ***********************************
                      .;;,
 .,.               .,;;;;;,
;;;;;;;,,        ,;;%%%%%;;
 `;;;%%%%;;,.  ,;;%%;;%%%;;
   `;%%;;%%%;;,;;%%%%%%%;;'
     `;;%%;;%:,;%%%%%;;%%;;,
        `;;%%%,;%%%%%%%%%;;;
           `;:%%%%%%;;%%;;;'
              .:::::::.
                   s.
 This routine performs the operation of a three dimensional complex to real
 inverse fft

     A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]

 Entry - A: a column distributed 3d block
         tmp2,tmp3: temporary work spaces that must be at
              least the size of (real) n2ft3d
 Exit - A is transformed and the imaginary
        part of A is set to zero
 uses - D3dB_c_transpose_jk, d3db_c_transpose_ijk
*/

void d3db_cr_fft3b(REAL A[], REAL tmp2[], REAL tmp3[])
{
   int i,j,k,q,indx, nxh,nxhy,nxhz,indx0,indx1;

   nxh = (nx/2+1);
   nxhz = nxh*nz;
   nxhy = nxh*ny;


   /**** slab mapping ****/
   if (mapping==1)
   {
      /* Do a transpose of A, A(kx,kz,ky) <- A(kx,ky,kz) */
      d3db_c_transpose_jk(A,tmp2,tmp3);


      /* fft along kz dimension, A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)] */
      indx0 = 0;
      for (q=0; q<nq; ++q)
      {
         for (i=0; i<nxh; ++i)
         {
            indx  = i + indx0;
            indx1 = indx;
            for (k=0; k<nz; ++k)
            {
               tmp2[2*k]   = A[2*indx];
               tmp2[2*k+1] = A[2*indx+1];
        	   indx += nxh;
            }
            ecfftb(&nz,tmp2,tmpz);
            for (k=0; k<nz; ++k)
            {
               A[2*indx1]   = tmp2[2*k];
               A[2*indx1+1] = tmp2[2*k+1];
               indx1 += nxh;
            }
         }
         indx0 += nxhz;
      }

      /* transpose of A, A(kx,ky,nz) <- A(kx,nz,ky) */
      d3db_c_transpose_jk(A,tmp2,tmp3);

      /* fft along ky dimension, A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] */
      indx0 = 0;
      for (q=0; q<nq; ++q)
      {
         for (i=0; i<(nx/2+1); ++i)
         {
            indx  = i + indx0;
            indx1 = indx;
            for (j=0; j<ny; ++j)
            {
               tmp2[2*j]   = A[2*indx];
               tmp2[2*j+1] = A[2*indx+1];
               indx += nxh;
            }
            ecfftb(&ny,tmp2,tmpy);
            for (j=0; j<ny; ++j)
            {
               A[2*indx1]   = tmp2[2*j];
               A[2*indx1+1] = tmp2[2*j+1];
               indx1 += nxh;
            }
         }
         indx0 += nxhy;
      }


      /* fft along kx dimension, A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] */
      cshift1_fftb(nx,ny,nq,1,A);
      indx = 0;
      for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
      {
         erfftb(&nx,&A[2*indx],tmpx);
         indx += nxh;
      }
      zeroend_fftb(nx,ny,nq,1,A);
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
       /* fft along kz dimension, A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] */
      d3db_fftbz_sub2(nq3,nz,tmpz,A);
      d3db_c_transpose_ijk(3,A,tmp2,tmp3);

      /* fft along ky dimension,A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] */
      d3db_fftby_sub2(nq2,ny,tmpy,A);
      d3db_c_transpose_ijk(4,A,tmp2,tmp3);

      /* fft along kx dimension, A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] */
      cshift1_fftb(nx,nq1,1,1,A);
      d3db_fftbx_sub(nq1,nx,nxh,tmpx,A);
      zeroend_fftb(nx,nq1,1,1,A);
   }
}



/***********************************
 *                                 *
 *         d3db_rc_fft3f           *
 *                                 *
 ***********************************
           ,,_
       zd$$??=
     z$$P? F:`c,                _
    d$$, `c'cc$$i           ,cd$?R
   $$$$ cud$,?$$$i       ,=P"2?z "
    $" " ?$$$,?$$$.    ,-''`>, bzP
     'cLdb,?$$,?$$$   ,h' "I$'J$P
  ... `?$$$,"$$,`$$h  $$PxrF'd$"
d$PP""?-,"?$$,?$h`$$,,$$'$F44"
?,,_`=4c,?=,"?hu?$`?L4$'? '
   `""?==""=-"" `""-`'_,,,,
           .ccu?m?e?JC,-,"=?
                     """=='?"
This routine performs the operation of a three dimensional real to complex
 fft

     A(kx,ky,kz) <- FFT3[A(nx,ny,nz)]

 Entry - A: a column distributed 3d block
         tmp2,tmp3: temporary work spaces that must be at
              least the size of (real) n2ft3d
 Exit - A is transformed and the padding is set to zero
 uses - d3dB_c_transpose_jk, d3db_c_transpose_ijk
*/

void d3db_rc_fft3f(REAL A[], REAL tmp2[], REAL tmp3[])
{
   int i,j,k,q,indx, nxh,nxhy,nxhz,indx0,indx1;

   nxh = (nx/2+1);
   nxhz = nxh*nz;
   nxhy = nxh*ny;

   /* slab mapping */
   if (mapping==1)
   {
      /* fft along nx dimension -  A(kx,ny,nz) <- fft1d[A(nx,ny,nz)] */
      indx = 0;
      for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
      {
         erfftf(&nx,&A[2*indx],tmpx);
    	 indx += nxh;
      }
      cshift_fftf(nx,ny,nq,1,A);


      /* fft along ny dimension -A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  */
      for (i=0; i<nxh; ++i)
      {
         indx  = i;
         indx1 = i;
         for (q=0; q<nq; ++q)
         {
            for (j=0; j<ny; ++j)
            {
        	   tmp2[2*j]   = A[2*indx];
        	   tmp2[2*j+1] = A[2*indx+1];
        	   indx += nxh;
            }

            ecfftf(&ny,tmp2,tmpy);

            for (j=0; j<ny; ++j)
            {
            	A[2*indx1]   = tmp2[2*j];
            	A[2*indx1+1] = tmp2[2*j+1];
            	indx1 += nxh;
            }
         }
      }

      /* transpose of A - A(ky,nz,ky) <- A(kx,ky,nz) */
      d3db_c_transpose_jk(A,tmp2,tmp3);

      /* fft along nz dimension - A(kx,kz,ky) <- fft1d[A(kx,nz,ky)] */
      for (i=0; i<nxh; ++i)
      {
    	  indx  = i;
    	  indx1 = i;
    	  for (q=0; q<nq; ++q)
    	  {
             for (k=0; k<nz; ++k)
             {
            	 tmp2[2*k]   = A[2*indx];
            	 tmp2[2*k+1] = A[2*indx+1];
            	 indx += nxh;
             }

             ecfftf(&nz,tmp2,tmpz);

             for (k=0; k<nz; ++k)
             {
                A[2*indx1]   = tmp2[2*k];
            	A[2*indx1+1] = tmp2[2*k+1];
            	indx1 += nxh;
             }
    	  }
      }

      /* transpose of A - A(kx,ky,kz) <- A(kx,kz,ky) */
     d3db_c_transpose_jk(A,tmp2,tmp3);
   }

   /* hilbert mapping */
   else
   {
      /* fft along nx dimension - A(kx,ny,nz) <- fft1d[A(nx,ny,nz)] */
	  indx = 0;
	  for (q=0; q<nq1; ++q)
	  {
		  erfftf(&nx,&A[2*indx],tmpx);
		  indx += nxh;
	  }
	  cshift_fftf(nx,nq1,1,1,A);

	  d3db_c_transpose_ijk(1,A,tmp2,tmp3);

	  /* fft along ny dimension - A(ky,nz,kx) <- fft1d[A(ny,nz,kx)] */
      indx = 0;
      for (q=0; q<nq2; ++q)
      {
    	  ecfftf(&ny,&A[2*indx],tmpy);
    	  indx += ny;
      }

      d3db_c_transpose_ijk(2,A,tmp2,tmp3);

      /* fft along nz dimension - A(kz,kx,ky) <- fft1d[A(nz,kx,ky)] */
      indx = 0;
      for (q=0; q<nq3; ++q)
      {
    	  ecfftf(&nz,&A[2*indx],tmpz);
    	  indx += nz;
      }

   }

}



/***********************************
 *					               *
 *	         d3db_cc_dot  	       *
 *					               *
 ***********************************

*/
REAL d3db_cc_dot(const REAL A[], const REAL B[])
{
    int i,j,k,q,index,p;
    int np = Parallel_np();
    int taskid = Parallel_taskid();
    REAL sum = 0.0;

   /**** slab mapping ****/
   if (mapping==1)
   {
      /* kx!=0 plane, so double count */
	  for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
	  for (i=1; i<(nx/2+1); ++i)
	  {
		  index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
		  sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  for (q=0; q<nq; ++q)
	  for (j=0; j<ny; ++j)
	  {
		 index = q*(nx/2+1)*ny + j*(nx/2+1);
		 sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
   }

   /* hilbert mapping */
   else
   {
      /* kx!=0 plane, so double count */
	  for (index=0; index<nfft3d_map; ++index)
	  {
		  sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  i = 0;
	  for (k=0; k<nz; ++k)
	  for (j=0; j<ny; ++j)
	  {
		  i = 0;
		  d3db_ijktoindexp(i,j,k,&index,&p);
		  if (p==taskid)
			  sum -= (A[2*index]*B[2*index]+A[2*index+1]*B[2*index+1]);
	  }
   }

   /* add up sums from other nodes */
   if (np>1)
      sum = Parallel_SumAll(sum);

   return sum;
}

/***********************************
 *					               *
 *	         d3db_cc_idot  	       *
 *					               *
 ***********************************

*/
REAL d3db_cc_idot(const REAL A[], const REAL B[])
{
    int i,j,k,q,index,p;
    int np = Parallel_np();
    int taskid = Parallel_taskid();
    REAL sum = 0.0;

   /**** slab mapping ****/
   if (mapping==1)
   {
      /* kx!=0 plane, so double count */
	  for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
	  for (i=1; i<(nx/2+1); ++i)
	  {
		  index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
		  sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  for (q=0; q<nq; ++q)
	  for (j=0; j<ny; ++j)
	  {
		 index = q*(nx/2+1)*ny + j*(nx/2+1);
		 sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
   }

   /* hilbert mapping */
   else
   {
      /* kx!=0 plane, so double count */
	  for (index=0; index<nfft3d_map; ++index)
	  {
		  sum += A[2*index]*B[2*index] + A[2*index+1]*B[2*index+1];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  i = 0;
	  for (k=0; k<nz; ++k)
	  for (j=0; j<ny; ++j)
	  {
		  i = 0;
		  d3db_ijktoindexp(i,j,k,&index,&p);
		  if (p==taskid)
			  sum -= (A[2*index]*B[2*index]+A[2*index+1]*B[2*index+1]);
	  }
   }

   return sum;
}


/***********************************
 *					               *
 *	         d3db_tt_dot  	       *
 *					               *
 ***********************************

*/
REAL d3db_tt_dot(const REAL A[], const REAL B[])
{
    int i,j,k,q,index,p;
    int np = Parallel_np();
    int taskid = Parallel_taskid();
    REAL sum = 0.0;

   /**** slab mapping ****/
   if (mapping==1)
   {
      /* kx!=0 plane, so double count */
	  for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
	  for (i=1; i<(nx/2+1); ++i)
	  {
		  index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
		  sum += A[index]*B[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  for (q=0; q<nq; ++q)
	  for (j=0; j<ny; ++j)
	  {
		 index = q*(nx/2+1)*ny + j*(nx/2+1);
		 sum += A[index]*B[index];
	  }
   }

   /* hilbert mapping */
   else
   {
      /* kx!=0 plane, so double count */
	  for (index=0; index<nfft3d_map; ++index)
	  {
		  sum += A[index]*B[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  i = 0;
	  for (k=0; k<nz; ++k)
	  for (j=0; j<ny; ++j)
	  {
		  i = 0;
		  d3db_ijktoindexp(i,j,k,&index,&p);
		  if (p==taskid)
			  sum -= (A[index]*B[index]);
	  }
   }

   /**** add up sums from other nodes ****/
   if (np>1)
	   sum = Parallel_SumAll(sum);

   return sum;
}



/***********************************
 *					               *
 *	         d3db_tt_idot  	       *
 *					               *
 ***********************************

*/
REAL d3db_tt_idot(const REAL A[], const REAL B[])
{
    int i,j,k,q,index,p;
    int np = Parallel_np();
    int taskid = Parallel_taskid();
    REAL sum = 0.0;

   /**** slab mapping ****/
   if (mapping==1)
   {
      /* kx!=0 plane, so double count */
	  for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
	  for (i=1; i<(nx/2+1); ++i)
	  {
		  index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
		  sum += A[index]*B[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  for (q=0; q<nq; ++q)
	  for (j=0; j<ny; ++j)
	  {
		 index = q*(nx/2+1)*ny + j*(nx/2+1);
		 sum += A[index]*B[index];
	  }
   }

   /* hilbert mapping */
   else
   {
      /* kx!=0 plane, so double count */
	  for (index=0; index<nfft3d_map; ++index)
	  {
		  sum += A[index]*B[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  i = 0;
	  for (k=0; k<nz; ++k)
	  for (j=0; j<ny; ++j)
	  {
		  i = 0;
		  d3db_ijktoindexp(i,j,k,&index,&p);
		  if (p==taskid)
			  sum -= (A[index]*B[index]);
	  }
   }

   return sum;
}


/***********************************
 *					               *
 *	         d3db_t_sum  	       *
 *					               *
 ***********************************

*/
REAL d3db_t_sum(const REAL A[])
{
    int i,j,k,q,index,p;
    int np = Parallel_np();
    int taskid = Parallel_taskid();
    REAL sum = 0.0;

   /**** slab mapping ****/
   if (mapping==1)
   {
      /* kx!=0 plane, so double count */
	  for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
	  for (i=1; i<(nx/2+1); ++i)
	  {
		  index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
		  sum += A[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  for (q=0; q<nq; ++q)
	  for (j=0; j<ny; ++j)
	  {
		 index = q*(nx/2+1)*ny + j*(nx/2+1);
		 sum += A[index];
	  }
   }

   /* hilbert mapping */
   else
   {
      /* kx!=0 plane, so double count */
	  for (index=0; index<nfft3d_map; ++index)
	  {
		  sum += A[index];
	  }
	  sum *= 2;

      /* kx==0 plane, so single count */
	  i = 0;
	  for (k=0; k<nz; ++k)
	  for (j=0; j<ny; ++j)
	  {
		  i = 0;
		  d3db_ijktoindexp(i,j,k,&index,&p);
		  if (p==taskid)
			  sum -= (A[index]);
	  }
   }

   /**** add up sums from other nodes ****/
   if (np>1)
	   sum = Parallel_SumAll(sum);

   return sum;
}




/***********************************
 *					               *
 *	         d3db_c_read  	       *
 *					               *
 ***********************************

*/
void d3db_c_read(FILE *fp, REAL A[], REAL tmp[], REAL tmp2[])
{
   int j,k,q,p_to,p_here;
   int taskid = Parallel_taskid();
   int np     = Parallel_np();
   int one = 1;
   int nn = (nx+2)*ny;
   int n  = (nx+2);
   int msgtype = 89;

   if (mapping==1)
   {
      /**** master node reads from file and distributes */
      if (taskid==0)
      {
         for (k=0; k<nz; ++k)
         {
        	 fread(tmp,sizeof(REAL),nn,fp);
        	 d3db_ktoqp(k,&q,&p_to);
             if (p_to==0)
                ecopy(&nn,tmp,&one,&A[q*nn],&one);
             else
                if (MPI_Send(tmp,nn,MPI_REAL_PRECISION,
            			      p_to,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
            		 printf("d3db_c_read error: MPI_Send failed\n");
         }
      }
      else
      {
    	  for (k=0; k<nz; ++k)
    	  {
    	     d3db_ktoqp(k,&q,&p_here);
    	     if (p_here==taskid)
                if (MPI_Recv(&A[q*nn],nn,MPI_REAL_PRECISION,0,msgtype,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_read error: MPI_Recv failed\n");
    	  }
      }
   }

   /* hilbert mapping */
   else
   {
      if (taskid==0)
      {
	     for (k=0; k<nz; ++k)
         for (j=0; j<ny; ++j)
         {
            fread(tmp,sizeof(REAL),n,fp);
            q    = q_map1[j+k*ny];
            p_to = p_map1[j+k*ny];
            if (p_to==0)
               ecopy(&n,tmp,&one,&A[q*n],&one);
            else
               if (MPI_Send(tmp,n,MPI_REAL_PRECISION,
            		     p_to,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
                  printf("d3db_c_read error: MPI_Send failed\n");
         }
      }
      else
      {
    	 for (k=0; k<nz; ++k)
         for (j=0; j<ny; ++j)
         {
             q      = q_map1[j+k*ny];
             p_here = p_map1[j+k*ny];
             if (p_here==taskid)
                if (MPI_Recv(&A[q*n],n,MPI_REAL_PRECISION,0,msgtype,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_read error: MPI_Recv failed\n");
      	  }
      }
      d3db_c_transpose_ijk(5,A,tmp,tmp2);
   }

}



/***********************************
 *					               *
 *	         d3db_c_write  	       *
 *					               *
 ***********************************

*/
void d3db_c_write(FILE *fp, REAL A[], REAL tmp[], REAL tmp2[])
{
   int j,k,q,p_here,p_from,idum;
   int taskid = Parallel_taskid();
   int np     = Parallel_np();
   int one = 1;
   int nn = (nx+2)*ny;
   int n  = (nx+2);
   int msgtype     = 89;
   int msgtype_dum = 87;

   if (mapping==1)
   {
      /**** master node write to file */
      if (taskid==0)
      {
         for (k=0; k<nz; ++k)
         {
        	 d3db_ktoqp(k,&q,&p_from);
             if (p_from==0)
                ecopy(&nn,&A[q*nn],&one,tmp,&one);
             else
             {
                if (MPI_Send(&idum,1,MPI_INTEGER,
            			      p_from,msgtype_dum,MPI_COMM_WORLD)!=MPI_SUCCESS)
            		 printf("d3db_c_write error: MPI_Send failed\n");
                if (MPI_Recv(tmp,nn,MPI_REAL_PRECISION,p_from,msgtype,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_read error: MPI_Recv failed\n");
             }
             fwrite(tmp,sizeof(REAL),nn,fp);
         }
      }
      else
      {
    	  for (k=0; k<nz; ++k)
    	  {
    	     d3db_ktoqp(k,&q,&p_here);
    	     if (p_here==taskid)
    	     {
                if (MPI_Recv(&idum,1,MPI_INTEGER,0,msgtype_dum,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_write error: MPI_Recv failed\n");
                if (MPI_Send(&A[q*nn],nn,MPI_REAL_PRECISION,0,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
                   printf("d3db_c_write error: MPI_Send failed\n");
    	     }
    	  }
      }
   }

   /* hilbert mapping */
   else
   {
      d3db_c_transpose_ijk(6,A,tmp,tmp2);

      if (taskid==0)
      {
    	 for (k=0; k<nz; ++k)
         for (j=0; j<ny; ++j)
         {
            q      = q_map1[j+k*ny];
            p_from = p_map1[j+k*ny];
            if (p_from==0)
            	ecopy(&n,&A[q*n],&one,tmp,&one);
            else
            {
                if (MPI_Send(&idum,1,MPI_INTEGER,
            			      p_from,msgtype_dum,MPI_COMM_WORLD)!=MPI_SUCCESS)
            		 printf("d3db_c_write error: MPI_Send failed\n");
                if (MPI_Recv(tmp,n,MPI_REAL_PRECISION,p_from,msgtype,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_read error: MPI_Recv failed\n");
             }
             fwrite(tmp,sizeof(REAL),n,fp);
         }
      }

      else
      {
     	 for (k=0; k<nz; ++k)
         for (j=0; j<ny; ++j)
         {
            q      = q_map1[j+k*ny];
            p_here = p_map1[j+k*ny];
            if (p_here==taskid)
            {
                if (MPI_Recv(&idum,1,MPI_INTEGER,0,msgtype_dum,MPI_COMM_WORLD,MPI_STATUS_IGNORE)!=MPI_SUCCESS)
                   printf("d3db_c_write error: MPI_Recv failed\n");
                if (MPI_Send(&A[q*n],n,MPI_REAL_PRECISION,0,msgtype,MPI_COMM_WORLD)!=MPI_SUCCESS)
                   printf("d3db_c_write error: MPI_Send failed\n");
            }
         }
      }
   }
}

