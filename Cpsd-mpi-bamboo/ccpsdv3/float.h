#ifndef _FLOAT_H_
#define _FLOAT_H_

/* single precision floats */
/*
#define REALTYPE  "single precision"
#define FMT1    "%f"
#define FMT8p3    "%8.3f"
#define FMTE8p3    "%8.3e"
#define FMT10p1    "%10.1f"
#define FMT12p3    "%12.3f"
#define FMTE12p3    "%12.3e"
#define FMTE13p5    "%13.5e"
#define FMTE15p5    "%15.5e"
#define FMTE18p7   "%18.7e"
#define FMTE19p10   "%19.10e"
#define FMTE20p10   "%20.10e"
#define FMT10   "%10.3f %10.3f %10.3f"
#define REAL    float
#define edot    sdot
#define escal   sscal
#define ecopy   scopy
#define eaxpy   saxpy
#define esyev   ssyev
#define egemm   sgemm
#define ieamax	isamax
*/

/* double precision floats */
#define REALTYPE  "double precision"
#define FMT1    "%lf"
#define FMT8p3    "%8.3lf"
#define FMTE8p3    "%8.3le"
#define FMT10p1    "%10.1lf"
#define FMT12p3    "%12.3lf"
#define FMTE12p3    "%12.3le"
#define FMTE13p5    "%13.5le"
#define FMTE15p5    "%15.5le"
#define FMTE18p7   "%18.7le"
#define FMTE19p10   "%19.10le"
#define FMTE20p10   "%20.10le"
#define FMT10   "%10.3lf %10.3lf %10.3lf"
#define REAL    double
#define edot    ddot_
#define escal   dscal_
#define ecopy   dcopy_
#define eaxpy   daxpy_
#define esyev   dsyev_
#define egemm   dgemm_
#define ieamax	idamax_


/* define BLAS routines */
extern REAL edot();
extern void escal();
extern void ecopy();

#endif
