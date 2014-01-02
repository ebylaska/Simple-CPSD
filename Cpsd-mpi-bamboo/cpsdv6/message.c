
/*     ================================================
*     THIS ROUTINE WRITES MESSAGE, TIME, AND DATE ON
*     I/O UNIT=10.  THIS ROUTINE IS MACHINE-DEPENDENT.
*     ================================================
*/
#include	<time.h>
#include	<stdio.h>

void message(const int n)
{
   int i;
   char *msg,date[26];
   time_t t0,t1;
   
   t1 = time(&t0);
   msg = ctime(&t0);
   for (i=0; i<25; ++i) 
      date[i] = msg[i];
   date[24] = '\0';
   date[25] = '\0';

   if (n==1) msg = "job started";
   if (n==2) msg = "iteration started";
   if (n==3) msg = "iteration ended";
   if (n==4) msg = "job completed";
   if (n==5) msg = "job terminated";

   printf("          >>>   %s on %s    <<<\n",msg,date);
}
