/*
 * get_g.c
 *
 *  Created on: Dec 5, 2013
 *      Author: bylaska
 *
 */

#include	"float.h"
#include	"Parallel.h"
#include	"d3db.h"

/**************************************************
 *                                                *
 *                 get_g                          *
 *                                                *
 **************************************************/

void get_g(REAL *unitg, REAL *g)
{
   int i,j,k,k1,k2,k3,nxh,nyh,nzh,indx,p;
   REAL gx,gy,gz;
   int taskid = Parallel_taskid();
   int nx     = d3db_nx();
   int ny     = d3db_ny();
   int nz     = d3db_nz();
   int nfft3d = d3db_nfft3d();
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;

   for (k3=(-nzh+1); k3<= nzh; ++k3)
   for (k2=(-nyh+1); k2<= nyh; ++k2)
   for (k1=0; k1<= nxh; ++k1)
   {
     gx = k1*unitg[0] + k2*unitg[3] + k3*unitg[6];
     gy = k1*unitg[1] + k2*unitg[4] + k3*unitg[7];
     gz = k1*unitg[2] + k2*unitg[5] + k3*unitg[8];
     i=k1; if (i<0) i += nx;
     j=k2; if (j<0) j += ny;
     k=k3; if (k<0) k += nz;
     d3db_ijktoindexp(i,j,k,&indx,&p);
     if (p==taskid)
     {
        //indx = i + j*(nffth+1) + k*(nffth+1)*nfft;
        g[indx]          = gx;
        g[indx+  nfft3d] = gy;
        g[indx+2*nfft3d] = gz;
     }
   }
}
