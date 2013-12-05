/* dirac_exchange.c - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the Dirac  
   exchange potential and energy of a spin density rho.

*/
#include	<stdio.h>
#include        <math.h>

#include	"float.h"
#include	"dirac_exchange.h"

static	REAL	alpha=0.6666666666666666667;


/* define constants */  
#define onethird   0.33333333333333333333
#define fourthird  1.33333333333333333333
#define twothird   0.66666666666666666667;
#define two_to_onethird 1.25992104989487319067
#define A               0.25992104989487319067
#define pi              3.14159265358979311600
#define fourpi         12.56637061435917246399

/********************************
 *			 	*
 *        set_Dirac_alpha	*
 *				*
 ********************************/

void set_Dirac_alpha(REAL aa)
{
   alpha = aa*twothird;
}

/********************************
 *				*
 *         Dirac_alpha		*
 *				*
 ********************************/

REAL Dirac_alpha()
{

   return alpha;
}


/********************************
 *				*
 *         Dirac_Exchange	*
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Dirac exchange functional.

   Entry - rho: the spin density
   Exit  - Vx:  the spin dependent exchange functional
	   exc:	  the energy density
*/

void  Dirac_Exchange(const int ispin, const int n2ft3d, const REAL rho[], REAL vx[], REAL ex[])
{
    int i;
    REAL n,xi;
    REAL n_onethird;
    REAL f,df;
    REAL ex_p, ex_f;
    REAL ux_p, ux_f;

   if (alpha<1.0e-9)
   {
     for (i=0; i<n2ft3d; ++i)
     {
        vx[i] = 0.0; 
        vx[i+(ispin-1)*n2ft3d] = 0.0;
        ex[i] = 0.0;
     }
   } 
   else if (ispin==1) 
   {
      for (i=0; i<n2ft3d; ++i)
      {
         n     = (rho[i] + rho[i]);
         n_onethird = pow((3.0*n/pi),onethird);
         xi = 0.0;
         ex[i] = -(9.0/8.0)*alpha*n_onethird;
         vx[i] = -(3.0/2.0)*alpha*n_onethird;
      }
   }
   else
   {
      for (i=0; i<n2ft3d; ++i)
      {
         n     = (rho[i] + rho[i+n2ft3d]);
         n_onethird = pow((3.0*n/pi),onethird);
   
         if (n > 0.0)
            xi = (rho[0] - rho[1])/n;
         else
            xi = 1.0;

         ex_p = -(9.0/8.0)*alpha*n_onethird;
         ux_p = -(3.0/2.0)*alpha*n_onethird;
     
         ex_f = two_to_onethird*ex_p;
         ux_f = fourthird*ex_f;
         f =  (  pow((1.0+xi),fourthird) 
               + pow((1.0-xi),fourthird) 
               - 2.0)/(2.0*A);

         df = fourthird*(  pow((1.0+xi),onethird)
                         - pow((1.0-xi),onethird))
                        /(2.0*A);

          ex[i] = ex_p + f*(ex_f - ex_p);

          vx[i] = ux_p + f*(ux_f-ux_p) + (+1.0-xi)*(df*(ex_f-ex_p));
          vx[i+n2ft3d] = ux_p + f*(ux_f-ux_p) + (-1.0-xi)*(df*(ex_f-ex_p));
       }
    }


} /* Dirac_Exchange */
