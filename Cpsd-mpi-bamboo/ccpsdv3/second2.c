#include	<time.h>
#include	"float.h"

REAL current_second()
{
  clock_t t0;
  REAL t;

  t0 = clock();
  t  = ((REAL) t0)/((REAL) CLOCKS_PER_SEC);
  return t;
}

