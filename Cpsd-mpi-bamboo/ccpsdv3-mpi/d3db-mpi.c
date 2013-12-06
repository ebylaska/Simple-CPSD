/*
 * d3db-mpi.c
 *
 *  Created on: Dec 5, 2013
 *      Author: bylaska
 */

#include <stdio.h>
#include <stdlib.h>
#include	"Parallel.h"

/********************/
/* static variables */
/********************/
static int nq,nx,ny,nz;
static int nfft3d,n2ft3d,nfft3d_map,n2ft3d_map;
static int zplane_size;
static int mapping, mapping2d;

/* slab mapping */
static int *q_map, *p_map, *k_map, *nqarray;
static int *i1_start, *i2_start, *iq_to_i1, *iq_to_i2;

/* hilbert/h-curve mapping */
static int *q_map1,*p_map1,nq1;
static int *q_map2,*p_map2,nq2;
static int *q_map3,*p_map3,nq3;
static int *h_i1_start[6],*h_i2_start[6],*h_iq_to_i1[6],*h_iq_to_i2[6];

/* fft workspace variables */
static double *tmpx,*tmpy,*tmpz;



/******************************************
 *                                        *
 *          generate_map_indexes          *
 *                                        *
 ******************************************/
static void generate_map_indexes(int taskid, int np, int ny, int nz,
		                         int p_map[], int q_map[], *nq_out)
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
		indx_q = q;
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
	   i = tmp_p[j+k*ny]
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
 *******************************************/
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
          ++p
          if (p>=np){p = 0; ++q;}
	   }
	   for (k=0; k<nz; ++k)
	      if (p_map[k]==taskid)
		     k_map[q_map[k]] = k;
	   nfft3d = (nx/2+1)*ny*nq
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
		   hilbert2d_map(ny,nz,    p_map1);
		   hilbert2d_map(nz,nx/2+1,p_map2);
		   hilbert2d_map(nx/2+1,ny,p_map3);
		}
		else
		{
			hcurve_map(ny,nz,    p_map1);
			hcurve_map(nz,nx/2+1,p_map2);
			hcurve_map(nx/2+1,ny,p_map3);
		}
		generate_map_indexes(taskid,np,ny,nz,    p_map1,q_map1,&nq1);
		generate_map_indexes(taskid,np,nz,nx/2+1,p_map2,q_map2,&nq2);
		generate_map_indexes(taskid,np,nx/2+1,ny,p_map3,q_map3,&nq3);

		nfft3d = (nx/2+1)*nq1;
		if (ny*nq2)>nfft3d) nfft3d = ny*nq2;
		if (nz*nq3)>nfft3d) nfft3d = nz*nq3;
		n2ft3d = 2*nfft3d;

		nfft3d_map = nz*nq3;
		n2ft3d_map = (nx+2)*nq1;
	}
}



/***********************************
 *					               *
 *	          d3db_qtok   		   *
 *					               *
 ***********************************/
void d3db_qtok(int q, int *k)
{
   *k = k_map[q];
}


/***********************************
 *					               *
 *	          d3db_ktoqp  		   *
 *					               *
 ***********************************/
void d3db_ktoqp(int k, int *q, int *p)
{
   *q = q_map[k];
   *p = p_map[k];
}

/***********************************
 *					               *
 *	       d3db_ijktoindexp	       *
 *					               *
 ***********************************/
void d3db_ijktoindexp(int i, int j, int k,
		              int *indx, int *p)
{
   int q;

   /**** slab mapping ***/
   if (mapping==1)
   {
      q = q_map[k];
      *p = p_map[k];
      *indx = i + j*(nx/2+1) + q*(nx/2+1)*ny;
   }

   /**** Hilbert mapping ****/
   else
   {
      q = q_map3[i+j*(nx/2+1)];
      *p = p_map3[i+j*(nx/2+1)];
      *indx = k + q*nz;
   }
}


/***********************************
 *                                 *
 *           d3db_ijktoindex1p     *
 *                                 *
 ***********************************/
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
 *           D3dB_ijktoindex2p     *
 *                                 *
 ***********************************/
void d3db_ijktoindex2p(int i,int j,int k,int *indx,int *p)
{
   int q;

   /**** slab mapping ****/
   if (mapping==1)
   {
      q = q_map[j];
      *p = p_map[j];
      *indx = i + k*(nx+2) + q*(nx+2)*ny;
   }

   /**** Hilbert mapping ****/
   else
   {
      q = q_map1[j+k*ny];
      *p = p_map1[j+k*ny];
      *indx = i + q*(nx+2);
   }
}


/***********************************
 *					               *
 *	        d3db_nfft3d		       *
 *					               *
 ***********************************/
int d3db_nfft3d()
{
   return nfft3d;
}


/***********************************
 *                                 *
 *         d3db_nfft3d_map         *
 *                                 *
 ***********************************/
int d3db_nfft3d_map()
{
   return nfft3d_map;
}


/***********************************
 *					               *
 *	        d3db_n2ft3d		       *
 *					               *
 ***********************************/
int d3db_n2ft3d()
{
   return n2ft3d;
}


/***********************************
 *                                 *
 *         d3db_n2ft3d_map         *
 *                                 *
 ***********************************/
int d3db_n2ft3d_map()
{
   return n2ft3d_map;
}



/***********************************
 *					               *
 *	        d3db_nq			       *
 *				            	   *
 ***********************************/
int d3db_nq()
{
   return nq;
}


/***********************************
 *					               *
 *	        d3db_nx			       *
 *					               *
 ***********************************/
int d3db_nx()
{
   return nx;
}


/***********************************
 *					               *
 *	        d3db_ny			       *
 *					               *
 ***********************************/

int d3db_ny()
{
   return ny;
}


/***********************************
 *					               *
 *	        d3db_nz			       *
 *					               *
 ***********************************/
int d3db_nz()
{
   return nz;
}

/************************************
 *                                 *
 *         d3db_zplane_size        *
 *                                 *
 ***********************************/
int d3db_zplane_size()
{
   return zplane_size;
}


/*****************************************
 *                                       *
 *        d3db_c_transpose_jk_init       *
 *                                       *
 *****************************************/
static void d3db_c_transpose_jk_init()
{
	int index1,index2,it,proc_to,proc_from,qhere,phere,qto,pto,qfrom,pfrom,itmp,i,j,k;

	int taskid = Parallel_taskid();
	int np = Parallel_np();

	iq_to_i1 = (int *) malloc((nx/2+1)*ny*nq,sizeof(int));
	iq_to_i2 = (int *) malloc((nx/2+1)*ny*nz*sizeof(int));
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
 *****************************************/
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
         pto   = p_map2[k+i*nz]
         qto   = q_map2[k+i*nz]

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
         pto   = p_map1[j+k*ny]
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
	tmpx = (double *) malloc(2*(2*nx+15)*sizeof(double));
	tmpy = (double *) malloc(2*(2*ny+15)*sizeof(double));
	tmpz = (double *) malloc(2*(2*nz+15)*sizeof(double));
	drffti(nx,tmpx);
	dcffti(ny,tmpy);
	dcffti(nz,tmpz);
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

/*****************************************
 *                                       *
 *              d3db_end                 *
 *                                       *
 *****************************************/

void d3dB_end()
{
	d3db_timereverse_end();
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
		for (int i=0; i<6; ++i)
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
 *****************************************/
void d3db_init(int nx_in, int ny_in, int nz_in, map_in)
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

	d3db_c_timereverse_init();
	d3db_fft_init();
}
