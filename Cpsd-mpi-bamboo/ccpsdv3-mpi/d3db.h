/*
 * d3db.h
 *
 *  Created on: Dec 9, 2013
 *      Author: bylaska
 */

#ifndef D3DB_H_
#define D3DB_H_

#include <stdio.h>
#include "float.h"

extern int d3db_nfft3d();
extern int d3db_nfft3d_map();
extern int d3db_n2ft3d();
extern int d3db_n2ft3d_map();
extern int d3db_nq();
extern int d3db_nx();
extern int d3db_ny();
extern int d3db_nz();
extern int d3db_zplane_size();

extern void d3db_end();
extern void d3db_init(int, int, int, int);
extern void d3db_cr_fft3b(REAL*, REAL*, REAL*);
extern void d3db_rc_fft3f(REAL*, REAL*, REAL*);
extern void d3db_r_zero_ends(REAL*);

REAL d3db_cc_dot(const REAL*, const REAL*);
REAL d3db_cc_idot(const REAL*, const REAL*);
REAL d3db_tt_dot(const REAL*, const REAL*);
REAL d3db_tt_idot(const REAL*, const REAL*);
REAL d3db_t_sum(const REAL*);

void d3db_c_read(FILE*, REAL*, REAL*, REAL*);
void d3db_t_read(FILE*, REAL*, REAL*, REAL*);
void d3db_c_write(FILE*, REAL*, REAL*, REAL*);

#endif /* D3DB_H_ */
