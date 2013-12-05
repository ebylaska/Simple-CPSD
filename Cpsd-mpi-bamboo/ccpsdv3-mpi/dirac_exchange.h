#ifndef	_DIRAC_EXHANGE_H_
#define _DIRAC_EXHANGE_H_

#include	"float.h"

extern void   set_Dirac_alpha(REAL aa);
extern REAL   Dirac_alpha();
extern void   Dirac_Exchange(const int ispin, const int n2ft3d, const REAL rho[], REAL vx[], REAL ex[]);


#endif
