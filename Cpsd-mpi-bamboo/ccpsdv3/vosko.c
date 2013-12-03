/* vosko.c
   Author - Eric Bylaska
*/

#include        <math.h>
#include	"float.h"
#include        "vosko.h"

static REAL beta=1.0;

/********************************
 *                              *
 *        set_Vosko_beta        *
 *                              *
 ********************************/
void set_Vosko_beta(REAL bb)
{
   beta = bb;
}



/*---- parameters given by vosko et al -----------------*/
#define ap  3.109070e-02 
#define af  1.554530e-02
#define x0p -1.049800e-01 
#define x0f -3.250000e-01
#define bp  3.727440e+00 
#define bf  7.060420e+00
#define cp  1.293520e+01 
#define cf  1.805780e+01
/*------------------------------------------------------*/

/*     constants calculated from vosko's parameters */
#define xp   -4.581653e-01  
#define xf   -5.772521e-01
#define qp    6.151991e+00  
#define qf    4.730927e+00
#define xx0p  1.255491e+01  
#define xx0f  1.586879e+01
#define cp1   3.109070e-02  
#define cf1   1.554530e-02
#define cp2   9.690228e-04  
#define cf2   2.247860e-03
#define cp3   1.049800e-01  
#define cf3   3.250000e-01
#define cp4   3.878329e-02  
#define cf4   5.249122e-02
#define cp5   3.075995e+00  
#define cf5   2.365463e+00
#define cp6   1.863720e+00  
#define cf6   3.530210e+00
#define dp1   6.218140e-02  
#define df1   3.109060e-02
#define dp2   1.938045e-03  
#define df2   4.495720e-03
#define dp3   1.049800e-01  
#define df3   3.250000e-01
#define dp4  -3.205972e-02  
#define df4  -1.779316e-02
#define dp5  -1.192972e-01  
#define df5  -1.241661e-01
#define dp6   1.863720e+00  
#define df6   3.530210e+00
#define dp7   9.461748e+00  
#define df7   5.595417e+00
#define fc    1.923661e+00  
#define fd    2.564881e+00
#define crs   7.876233e-01

#define	small_number	1.0e-30

/********************************
 *				*
 *   correlation_vosko et. al.  *
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Vosko et. al. correlation functional.
   This is a Ceperly and Alder parameterization

   Entry - ispin:
           n2ft3d:
           dn[]: the spin density
   Exit  - xcp[]:  the spin dependent exchange-correlation functional
	   xce[]:  the exchange-correlation energy density
*/
   /* define constants */  
#define onesixth	0.166666666666667
#define onethird	0.333333333333333
#define fourthird	1.333333333333333
#define twothird	0.666666666666667
#define two_to_onethird	1.259921049894873
#define A               0.259921049894873
#define pi		3.141592653589793
#define rs_scale	0.620350490899400

void  Vosko(const int ispin, const int n2ft3d, const REAL dn[], REAL xcp[], REAL *xce)
{
    int i;
    REAL rs,nup,ndown,n,xi;
    REAL x,xxp,dxxp,xxf,dxxf;
    REAL f,df;
    REAL denominator;
    REAL ec_p, ec_f;
    REAL uc_p, uc_f;


   /* define constants */  
/*
   onesixth  = 1.0/6.0;
   onethird  = 1.0/3.0;
   fourthird = 4.0/3.0;
   twothird  = 2.0/3.0;
   two_to_onethird = pow(2.0,onethird);
   A               = two_to_onethird - 1.0;
   pi       = 4.0*atan(1.0);
   rs_scale = pow( (0.75/pi), (1.0/3.0));
*/

   if (beta>1.0e-9)
   for (i=0; i<n2ft3d; ++i)
   {
      nup   = dn[i]                  + small_number;
      ndown = dn[i+(ispin-1)*n2ft3d] + small_number;

      n     = (nup + ndown);
      rs    = rs_scale/pow(n,onethird);
      xi    = (nup - ndown)/n;

      x     = sqrt(rs);

      xxp  = rs + bp*x + cp;
      dxxp = 2.0*x + bp;
      ec_p = cp1*log(rs/xxp) + cp2*log( (x+cp3)*(x+cp3)/xxp)
           + cp4*atan(cp5/(x+cp6));
      uc_p = ec_p
            - onesixth*x*(  dp1/x + dp2/(x+cp3) + dp4*dxxp/xxp
                          + dp5/( (x+dp6)*(x+dp6)+dp7) 
                         );

      xxf  = rs + bf*x + cf;
      dxxf = 2.0*x + bf;
      ec_f = cf1*log(rs/xxf) + cf2*log( (x+cf3)*(x+cf3)/xxf)
           + cf4*atan(cf5/(x+cf6));
      uc_f = ec_f
            - onesixth*x*(  df1/x + df2/(x+cf3) + df4*dxxf/xxf
                          + df5/( (x+df6)*(x+df6)+df7) 
                         );
               
      f =  (  pow((1.0+xi),fourthird) 
            + pow((1.0-xi),fourthird) 
            - 2.0)/(2.0*A);

      df = fourthird*(  pow((1.0+xi),onethird)
                      - pow((1.0-xi),onethird))
                     /(2.0*A);

       xce[i] += beta*(ec_p + f*(ec_f - ec_p));

       xcp[i] += beta*(uc_p + f*(uc_f-uc_p) + (+1.0-xi)*(df*(ec_f-ec_p)));
       if (ispin==2)
          xcp[i+n2ft3d] += beta*(uc_p + f*(uc_f-uc_p) + (-1.0-xi)*(df*(ec_f-ec_p)));
   }
    
} /* Vosko */
