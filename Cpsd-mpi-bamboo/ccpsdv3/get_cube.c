    
#include	<stdlib.h>
#include	<stdio.h>  
#include        <math.h>

#include	"float.h"

/*************************************************
 *                                               *
 *                get_cube                       *
 *                                               *
 *************************************************/

/*  This routine computes primitive vectors both in coordination      
      space and in reciporocal space and the volume of primitive cell. 
                                                              
      Input:  itype --- type of cube (1=SC, 2=FCC, 3=BCC)     
               unit --- lattice constant                        

      Output: volume --- volume of primitive cell                  
              unita  --- primitive vectors in coordination space    
              unitg  --- primitive vectors in reciprocal space       

      Library:  DSCAL from BLAS                                        

      Last modification:  7/03/93  by R. Kawai                           
*/

void get_cube(int itype, REAL unit, REAL *volume, REAL *unita, REAL *unitg)
{
   REAL s,v;
   int nine=9;
   int one=1;
   REAL twopi = 8.0*atan(1.0);
 

   /* primitive vectors in coordination space */

   /*simple cubic*/
   if (itype==1)
   {
      unita[0] = 1.00; unita[1] = 0.00; unita[2] = 0.00;
      unita[3] = 0.00; unita[4] = 1.00; unita[5] = 0.00;
      unita[6] = 0.00; unita[7] = 0.00; unita[8] = 1.00;

   }
   /* face centered cubic */
   else if (itype==2)
   {
      unita[0] = +0.50; unita[1] = +0.50; unita[2] = +0.00;
      unita[3] = +0.50; unita[4] = +0.00; unita[5] = +0.50;
      unita[6] = +0.00; unita[7] = +0.50; unita[8] = +0.50;
 
   }
   /* body centered cubic */
   else if (itype==3)
   {
      unita[0] = +0.50; unita[1] = +0.50; unita[2] = +0.50;
      unita[3] = -0.50; unita[4] = +0.50; unita[5] = +0.50;
      unita[6] = +0.50; unita[7] = +0.50; unita[8] = -0.50;
   }
   /* unknown lattice itype */
   else
     printf("get_cube: input error -- unknown lattice itype.\n");

   escal(&nine,&unit,unita,&one);

   /* primitive vectors in the reciprocal space */
   unitg[0+0*3] = unita[1+1*3]*unita[2+2*3] - unita[2+1*3]*unita[1+2*3];
   unitg[1+0*3] = unita[2+1*3]*unita[0+2*3] - unita[0+1*3]*unita[2+2*3];
   unitg[2+0*3] = unita[0+1*3]*unita[1+2*3] - unita[1+1*3]*unita[0+2*3];
   unitg[0+1*3] = unita[1+2*3]*unita[2+0*3] - unita[2+2*3]*unita[1+0*3];
   unitg[1+1*3] = unita[2+2*3]*unita[0+0*3] - unita[0+2*3]*unita[2+0*3];
   unitg[2+1*3] = unita[0+2*3]*unita[1+0*3] - unita[1+2*3]*unita[0+0*3];
   unitg[0+2*3] = unita[1+0*3]*unita[2+1*3] - unita[2+0*3]*unita[1+1*3];
   unitg[1+2*3] = unita[2+0*3]*unita[0+1*3] - unita[0+0*3]*unita[2+1*3];
   unitg[2+2*3] = unita[0+0*3]*unita[1+1*3] - unita[1+0*3]*unita[0+1*3];

   v = unita[0]*unitg[0] + unita[1]*unitg[1] + unita[2]*unitg[2];
   s = twopi/v;
   escal(&nine,&s,unitg,&one);


   *volume=fabs(v);
}

