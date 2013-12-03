#include	<sys/time.h>
#include	"float.h"

#define bigint ((long) 1236025000)

REAL current_second()
{
  struct timeval  tp;
  struct timezone tz;
  double t0,t1,t2;

  gettimeofday(&tp,&tz); 
  t0  = ((double) (tp.tv_sec-bigint));
  t0 += ((double) tp.tv_usec)*1e-6;

  return ((REAL) t0);
}

