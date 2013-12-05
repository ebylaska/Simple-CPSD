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




static void generate_map_indexes(int taskid, int np, int ny, int nz, int p_map[], int q_map[], *nq_out)
{
	int i,nq1,
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
}

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

static void d3db_c_transpose_jk_init()

static void d3db_c_transpose_ijk_init()

static void d3db_fft_init()
{
	tmpx = (double *) malloc(2*(2*nx+15)*sizeof(double));
	tmpy = (double *) malloc(2*(2*ny+15)*sizeof(double));
	tmpz = (double *) malloc(2*(2*nz+15)*sizeof(double));
	drffti(nx,tmpx);
	dcffti(ny,tmpy);
	dcffti(nz,tmpz);
}
static void d3db_fft_end()
{
	free(tmpx);
	free(tmpy);
	free(tmpz);
}

void d3dB_end()
{
	D3dB_timereverse_end();

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
