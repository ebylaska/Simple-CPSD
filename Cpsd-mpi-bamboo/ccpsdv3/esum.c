#include	"float.h"

REAL esum(int n, REAL *sx, int incx)
{     
   int  i;
   REAL stemp;

    stemp = 0.0;
    for (i=0; i<n; ++i)
      stemp += sx[i*incx];

   return (stemp);
}

