    
#include	<math.h>
#include	<stdio.h>  

#define	REAL	double



REAL simp(int n, REAL *Y, REAL H)
{
   int ne,no;
   REAL S;
   ne = n/2;
   no = ne+1;

   S= 2.0*esum(no,Y,2) 
    + 4.0*esum(ne,&Y[1],2)
    -Y[0]-Y[n-1];
   
  return (S*H/3.0);
}

REAL esum(int n, double *sx, int incx)
{
   int  i;
   REAL stemp;

    stemp = 0.0;
    for (i=0; i<n; i+=incx)
      stemp += sx[i];

   return (stemp);
}
